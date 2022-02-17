/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_closedforms.cpp,v $
 * Revision 1.1  2004/02/05 15:08:43  ebenhamou,ocroissant
 * Initial version
 *
 */

/*! \file ARM_xl_gp_closedforms.cpp
 *
 *  \brief file for the closed forms part in the generic pricer
 *
 *	\author  O Croissant
 *	\version 1.0
 *	\date February 2004
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_gp_closedforms.h>
#include "ARM_xl_gp_closedforms_local.h"
#include "ARM_xl_wrapper_local.h"
#include <ARM\libarm_local\ARM_local_wrapper.h>

#include "ARM_gp_local_interglob.h"

#include "zerocurv.h"
#include <deque>

#include "ARM_xl_gp_fctorhelper.h"

#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

////////////////////////////////////////////
/// Vanilla Option Normal
////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaOption_Normal(
	LPXLOPER XL_Undelying,
	LPXLOPER XL_Volatility,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_CallOrPut)
{
	ADD_LOG("Local_VanillaOption_Normal");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs
	    double C_underling;
	    XL_readNumCell(XL_Undelying,C_underling,		" ARM_ERR: undelying numeric expected",C_result);	
	    double C_volatility;
	    XL_readNumCell(XL_Volatility,C_volatility,		" ARM_ERR: volatility numeric expected",C_result);	
	    double C_strike;
	    XL_readNumCell(XL_Strike,C_strike,	            " ARM_ERR: strike numeric expected",C_result);	
	    double C_maturity;
	    XL_readNumCell(XL_Maturity,C_maturity,	        " ARM_ERR: maturity numeric expected",C_result);	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);

	    /// call the function
	    long retCode = ARMLOCAL_VanillaOption_Normal(
		    C_underling,
		    C_volatility,
		    C_strike,
		    C_maturity,
		    C_callOrPut,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaOption_Normal_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_Undelying,
	LPXLOPER XL_Volatility,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_CallOrPut)
{
	ADD_LOG("Local_VanillaOption_Normal_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs
		double C_i;
		XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);
	    double C_underling;
	    XL_readNumCell(XL_Undelying,C_underling,		" ARM_ERR: undelying numeric expected",C_result);	
	    double C_volatility;
	    XL_readNumCell(XL_Volatility,C_volatility,		" ARM_ERR: volatility numeric expected",C_result);	
	    double C_strike;
	    XL_readNumCell(XL_Strike,C_strike,	            " ARM_ERR: strike numeric expected",C_result);	
	    double C_maturity;
	    XL_readNumCell(XL_Maturity,C_maturity,	        " ARM_ERR: maturity numeric expected",C_result);	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);

	    /// call the function
	    long retCode = ARMLOCAL_VanillaOption_Normal_Der(
			C_i,
		    C_underling,
		    C_volatility,
		    C_strike,
		    C_maturity,
		    C_callOrPut,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_VanillaOption_Normal_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_Undelying,
	LPXLOPER XL_Volatility,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_CallOrPut)
{
	ADD_LOG("Local_VanillaOption_Normal_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_i;
		XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);
		double C_j;
		XL_readNumCell(XL_j,C_j,		" ARM_ERR: rank numeric expected",C_result);
		double C_underling;
		XL_readNumCell(XL_Undelying,C_underling,		" ARM_ERR: undelying numeric expected",C_result);	
	    double C_volatility;
	    XL_readNumCell(XL_Volatility,C_volatility,		" ARM_ERR: volatility numeric expected",C_result);	
	    double C_strike;
	    XL_readNumCell(XL_Strike,C_strike,	            " ARM_ERR: strike numeric expected",C_result);	
	    double C_maturity;
	    XL_readNumCell(XL_Maturity,C_maturity,	        " ARM_ERR: maturity numeric expected",C_result);	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);


	    /// call the function
	    long retCode = ARMLOCAL_VanillaOption_Normal_Der2(
			C_i,
			C_j,
		    C_underling,
		    C_volatility,
		    C_strike,
		    C_maturity,
		    C_callOrPut,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_DoubleDigital_Normal(
	LPXLOPER XL_Fwd1,
	LPXLOPER XL_Fwd2,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_Strike1,
	LPXLOPER XL_Spread1,
	LPXLOPER XL_Strike2,
	LPXLOPER XL_Spread2,
	LPXLOPER XL_Vol1plus,
	LPXLOPER XL_Vol1minus,
	LPXLOPER XL_Vol2plus,
	LPXLOPER XL_Vol2minus,
	LPXLOPER XL_Correl,
	LPXLOPER XL_CallOrPut1,
	LPXLOPER XL_CallOrPut2)
{
	ADD_LOG("Local_DoubleDigital_Normal");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs
	    double C_Fwd1;
		XL_readNumCell(XL_Fwd1, C_Fwd1,			 " ARM_ERR: fwd1 numeric expected",C_result);	
		double C_Fwd2;
		XL_readNumCell(XL_Fwd2, C_Fwd2,			 " ARM_ERR: fwd2 numeric expected",C_result);	
		double C_Maturity;
		XL_readNumCell(XL_Maturity, C_Maturity,	 " ARM_ERR: maturity numeric expected",C_result);	
		double C_Strike1;
		XL_readNumCell(XL_Strike1,  C_Strike1,	 " ARM_ERR: strike1 numeric expected",C_result);	
		double C_Spread1;
		XL_readNumCell(XL_Spread1,  C_Spread1,	 " ARM_ERR: spread1 numeric expected",C_result);	
		double C_Strike2;
		XL_readNumCell(XL_Strike2,  C_Strike2,	 " ARM_ERR: strike2 numeric expected",C_result);	
		double C_Spread2;
		XL_readNumCell(XL_Spread2,  C_Spread2,	 " ARM_ERR: spread2 numeric expected",C_result);	
		double C_Vol1plus;
		XL_readNumCell(XL_Vol1plus, C_Vol1plus,	 " ARM_ERR: vol1plus numeric expected",C_result);	
		double C_Vol1minus;
		XL_readNumCell(XL_Vol1minus,C_Vol1minus, " ARM_ERR: vol1minus numeric expected",C_result);	
		double C_Vol2plus;
		XL_readNumCell(XL_Vol2plus, C_Vol2plus,	 " ARM_ERR: vol2plus numeric expected",C_result);	
		double C_Vol2minus;
		XL_readNumCell(XL_Vol2minus,C_Vol2minus, " ARM_ERR: vol2minus numeric expected",C_result);	
		double C_Correl;
		XL_readNumCell(XL_Correl,   C_Correl,    " ARM_ERR: correl numeric expected",C_result);
		
		double C_CallOrPut1;
        XL_GETCONVCALLORPUT(XL_CallOrPut1, C_CallOrPut1,   " ARM_ERR: call or put 1 string expected",C_result);
		double C_CallOrPut2;
        XL_GETCONVCALLORPUT(XL_CallOrPut2, C_CallOrPut2,   " ARM_ERR: call or put 2 string expected",C_result);
		

	    /// call the function
	    long retCode = ARMLOCAL_DoubleDigital_Normal (
									C_Fwd1, 
									C_Fwd2, 
									C_Maturity, 
									C_Strike1, 
									C_Spread1, 
									C_Strike2, 
									C_Spread2, 
									C_Vol1plus, 
									C_Vol1minus, 
									C_Vol2plus, 
									C_Vol2minus, 
									C_Correl, 
									C_CallOrPut1, 
									C_CallOrPut2,
									C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DoubleDigital_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




////////////////////////////////////////////
/// Spread Option Log Normal
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
{
	ADD_LOG("Local_LogNormal_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_SpreadOption(
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
	{
		ADD_LOG("Local_Smiled_LogNormal_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Smiled_LogNormal_SpreadOption(
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_LogNormal_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/// Calibration


__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption_Calibrate_Correlation(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
{
	ADD_LOG("Local_LogNormal_SpreadOption_Calibrate_Correlation");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_SpreadOption_Calibrate_Correlation(
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_SpreadOption_Calibrate_Correlation" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
	{
		ADD_LOG("Local_Smiled_LogNormal_SpreadOption_Calibrate_Correlation");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_LogNormal_SpreadOption_Calibrate_Correlation" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
{
	ADD_LOG("Local_LogNormal_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
	{
		ADD_LOG("Local_Smiled_LogNormal_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();	

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Smiled_LogNormal_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_LogNormal_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
{
	ADD_LOG("Local_LogNormal_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);	
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n )
	{
		ADD_LOG("Local_Smiled_LogNormal_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD( XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Smiled_LogNormal_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_LogNormal_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////
/// Spread Option Normal
////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype)
{
	ADD_LOG("Local_Normal_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Normal_SpreadOption(
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_Normal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype)
	{
		ADD_LOG("Local_Smiled_Normal_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	
	/// call the function
	long retCode=ARMLOCAL_Smiled_Normal_SpreadOption(
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_Normal_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Normal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype)
{
	ADD_LOG("Local_Normal_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Normal_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_Normal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype)
	{
		ADD_LOG("Local_Smiled_Normal_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Smiled_Normal_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_Normal_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Normal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype)
{
	ADD_LOG("Local_Normal_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);	
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Normal_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_callput,
		C_optiontype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_Normal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype)
	{
		ADD_LOG("Local_Smiled_Normal_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: rank numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: rank numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_slope1;
	XL_readNumCell(XL_slope1,C_slope1,		" ARM_ERR: t numeric expected",C_result);
	double C_slope2;
	XL_readNumCell(XL_slope2,C_slope2,		" ARM_ERR: t numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optiontype;
	XL_GETCONVSPREADDERIVEDSTRUCT(XL_optiontype,C_optiontype," ARM_ERR: SpreadOption Derived Structure expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_Smiled_Normal_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_sig1,
		C_sig2,
		C_rho,
		C_k,
		C_t,
		C_slope1,
		C_slope2,
		C_callput,
		C_optiontype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Smiled_Normal_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// Spread Option Gaussian
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_Parameter_Corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_Parameter_Corr,C_copula_corr,	" ARM_ERR: correlation expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of 4 numeric expected",C_result);




	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_SpreadOption(
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_t,
		C_flag,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// Spread Option SABR Gaussian
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_Digital_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_Digital_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption(
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_t,
		C_flag,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_t,
		C_flag,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_Digital_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_Digital_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_t,
		C_flag,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_t,
		C_flag,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_Digital_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_Digital_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_t,
		C_flag,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption_Certitude(
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Gaussian_SABR_Power_SpreadOption_Certitude");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Certitude(
		C_copula_corr,
		C_t,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption_Certitude" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Formula
///
////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_CF_LN_VanillaOption(
	LPXLOPER XL_forward,
	LPXLOPER XL_volatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_maturity,
	LPXLOPER XL_CallPut
	)
	{
		ADD_LOG("Local_CF_LN_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_volatility;
	XL_readNumCell(XL_volatility,C_volatility,		" ARM_ERR: totalvolatility numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_maturity;
	XL_readNumCell(XL_maturity,C_maturity,	" ARM_ERR: maturity numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholes(
		C_forward,
		C_volatility*sqrt(C_maturity),
		C_bondprice,
		C_strike,
		C_callput,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_LN_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes(
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	)
	{
		ADD_LOG("Local_CF_BlackSholes");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_totalvolatility;
	XL_readNumCell(XL_totalvolatility,C_totalvolatility,		" ARM_ERR: totalvolatility numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholes(
		C_forward,
		C_totalvolatility,
		C_bondprice,
		C_strike,
		C_callput,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	)
	{
		ADD_LOG("Local_CF_BlackSholes_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_totalvolatility;
	XL_readNumCell(XL_totalvolatility,C_totalvolatility,		" ARM_ERR: totalvolatility numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholes_Der(
		C_i,
		C_forward,
		C_totalvolatility,
		C_bondprice,
		C_strike,
		C_callput,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	)
	{
		ADD_LOG("Local_CF_BlackSholes_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_totalvolatility;
	XL_readNumCell(XL_totalvolatility,C_totalvolatility,		" ARM_ERR: totalvolatility numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholes_Der2(
		C_i,
		C_j,
		C_forward,
		C_totalvolatility,
		C_bondprice,
		C_strike,
		C_callput,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_ImplicitVolatility(
	LPXLOPER XL_forward,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_maturity,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_optprice,
	LPXLOPER XL_algo
	)
	{
		ADD_LOG("Local_CF_BlackSholes_ImplicitVolatility");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: strike numeric expected",C_result);	
	double C_maturity;
	XL_readNumCell(XL_maturity,C_maturity,	" ARM_ERR: maturity numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optprice;
	XL_readNumCell(XL_optprice,C_optprice,		" ARM_ERR: forward numeric expected",C_result);	
	double C_algo;
	double C_algo_Default=1;
	XL_readNumCellWD(XL_algo,C_algo,C_algo_Default,		" ARM_ERR: algo numeric expected",C_result);	


	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholes_ImplicitVolatility(
		C_forward,
		C_bondprice,
		C_strike,
		C_maturity,
		C_callput,
		C_optprice,
		C_algo,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes_ImplicitVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_ImplicitVol(
	LPXLOPER XL_forward,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_optprice,
	LPXLOPER XL_algo
	)
	{
		ADD_LOG("Local_CF_BlackSholes_ImplicitVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optprice;
	XL_readNumCell(XL_optprice,C_optprice,		" ARM_ERR: forward numeric expected",C_result);	
	double C_algo;
	double C_algo_Default=1;
	XL_readNumCellWD(XL_algo,C_algo,C_algo_Default,		" ARM_ERR: algo numeric expected",C_result);	


	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholes_ImplicitVol(
		C_forward,
		C_bondprice,
		C_strike,
		C_callput,
		C_optprice,
		C_algo,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes_ImplicitVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
/// Spread Option ShiftedLognormal Gaussian
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Gaussian_ShiftedLN_Power_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sigma1;
	XL_readNumCell(XL_sigma1,C_sigma1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_sigma2;
	XL_readNumCell(XL_sigma2,C_sigma2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption(
		C_S1,
		C_S2,
		C_sigma1,
		C_alpha1,
		C_sigma2,
		C_alpha2,
		C_copula_corr,
		C_t,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_ShiftedLN_Power_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Gaussian_ShiftedLN_Power_SpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sigma1;
	XL_readNumCell(XL_sigma1,C_sigma1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_sigma2;
	XL_readNumCell(XL_sigma2,C_sigma2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_sigma1,
		C_alpha1,
		C_sigma2,
		C_alpha2,
		C_copula_corr,
		C_t,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_ShiftedLN_Power_SpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Gaussian_ShiftedLN_Power_SpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_sigma1;
	XL_readNumCell(XL_sigma1,C_sigma1,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_sigma2;
	XL_readNumCell(XL_sigma2,C_sigma2,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_sigma1,
		C_alpha1,
		C_sigma2,
		C_alpha2,
		C_copula_corr,
		C_t,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_ShiftedLN_Power_SpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Gaussian_ShiftedLN_Power_SpreadOption_Certitude");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
		C_copula_corr,
		C_t,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_ShiftedLN_Power_SpreadOption_Certitude" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// Merton Jump Diffusion Formula
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Merton_JumpDiffusion(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_n
	)
	{
		ADD_LOG("Local_Merton_JumpDiffusion");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);	
	double C_sigma;
	XL_readNumCell(XL_sigma,C_sigma,		" ARM_ERR: sigma numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	double C_muJ;
	XL_readNumCell(XL_muJ,C_muJ,	" ARM_ERR: muJ numeric expected",C_result);	
	double C_sigmaJ;
	XL_readNumCell(XL_sigmaJ,C_sigmaJ,	" ARM_ERR: sigmaJ numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Merton_JumpDiffusion(
		C_F,
		C_K,
		C_t,
		C_sigma,
		C_lambda,
		C_muJ,
		C_sigmaJ,
		C_callput,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Merton_JumpDiffusion" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Merton_JumpDiffusion_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_n
	)
	{
		ADD_LOG("Local_Merton_JumpDiffusion_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);	
	double C_sigma;
	XL_readNumCell(XL_sigma,C_sigma,		" ARM_ERR: sigma numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	double C_muJ;
	XL_readNumCell(XL_muJ,C_muJ,	" ARM_ERR: muJ numeric expected",C_result);	
	double C_sigmaJ;
	XL_readNumCell(XL_sigmaJ,C_sigmaJ,	" ARM_ERR: sigmaJ numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Merton_JumpDiffusion_Der(
		C_i,
		C_F,
		C_K,
		C_t,
		C_sigma,
		C_lambda,
		C_muJ,
		C_sigmaJ,
		C_callput,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Merton_JumpDiffusion_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Merton_JumpDiffusion_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_n
	)
	{
		ADD_LOG("Local_Merton_JumpDiffusion_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);	
	double C_sigma;
	XL_readNumCell(XL_sigma,C_sigma,		" ARM_ERR: sigma numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	double C_muJ;
	XL_readNumCell(XL_muJ,C_muJ,	" ARM_ERR: muJ numeric expected",C_result);	
	double C_sigmaJ;
	XL_readNumCell(XL_sigmaJ,C_sigmaJ,	" ARM_ERR: sigmaJ numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Merton_JumpDiffusion_Der2(
		C_i,
		C_j,
		C_F,
		C_K,
		C_t,
		C_sigma,
		C_lambda,
		C_muJ,
		C_sigmaJ,
		C_callput,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Merton_JumpDiffusion_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR Implicit Vol Formula
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_SABR_ImplicitVol(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_ImplicitVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_SABR_ImplicitVol(
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_ImplicitVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_sabr2b_ImplicitVol(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta1,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_zero,
	LPXLOPER XL_lambda)
{
	ADD_LOG("Local_sabr2b_ImplicitVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_zero;
	XL_readNumCell(XL_zero,C_zero,	" ARM_ERR: zero numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	
	/// call the function
	long retCode=ARMLOCAL_SABR2B_ImplicitVol(
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta1,
		C_beta2,
		C_rho,
		C_nu,
		C_zero,
		C_lambda,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_ImplicitVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_ImplicitVol_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_ImplicitVol_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_ImplicitVol_Der(
		C_i,
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_ImplicitVol_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_ImplicitVol_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_ImplicitVol_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
		double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_ImplicitVol_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_ImplicitVol_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SABR_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_callput,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_VanillaOption(
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_callput,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_FromGaussianToDistribution(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_FromGaussianToDistribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_SABR_FromGaussianToDistribution(
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_SABR_FromDistributionToGaussian(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_FromDistributionToGaussian");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_SABR_FromDistributionToGaussian(
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_callput,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_VanillaOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_VanillaOption_Der(
		C_i,
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_callput,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_VanillaOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_callput,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_VanillaOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_VanillaOption_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_tex,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_callput,
		C_flag,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_VanillaOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_From_Sigma_To_Alpha(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_sigma,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps
	)
	{
		ADD_LOG("Local_SABR_From_Sigma_To_Alpha");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_tex;
	XL_readNumCell(XL_tex,C_tex,		" ARM_ERR: tex numeric expected",C_result);	
	double C_sigma;
	XL_readNumCell(XL_sigma,C_sigma,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,	" ARM_ERR: nu numeric expected",C_result);	
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_nbsteps;
	double n_default;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_From_Sigma_To_Alpha(
		C_f,
		C_K,
		C_tex,
		C_sigma,
		C_beta,
		C_rho,
		C_nu,
		C_flag,
		C_nbsteps,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_From_Sigma_To_Alpha" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///		Single Barriere Option in Black and Sholes Model Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_discount,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	)
	{
		ADD_LOG("Local_BS_EuroBarriere");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_discount;
	XL_readNumCell(XL_discount,C_discount,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_inout;
	XL_GETCONVBARRIEREINOUT(XL_inout,C_inout," ARM_ERR: IN or OUT flag string expected",C_result);
	double C_updown;
	XL_GETCONVBARRIEREUPDOWN(XL_updown,C_updown," ARM_ERR: UP or DOWN flag string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroBarriere(
		C_f,
		C_k,
		C_b,
		C_r,
		C_v,
		C_t,
		C_discount,
		C_callput,
		C_inout,
		C_updown,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroBarriere" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere_ImpliedVol(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_op,
	LPXLOPER XL_t,
	LPXLOPER XL_discount,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	)
	{
		ADD_LOG("Local_BS_EuroBarriere_ImpliedVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_op;
	XL_readNumCell(XL_op,C_op,	" ARM_ERR: op numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_discount;
	XL_readNumCell(XL_discount,C_discount,	" ARM_ERR: discount numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_inout;
	XL_GETCONVBARRIEREINOUT(XL_inout,C_inout," ARM_ERR: IN or OUT flag string expected",C_result);
	double C_updown;
	XL_GETCONVBARRIEREUPDOWN(XL_updown,C_updown," ARM_ERR: UP or DOWN flag string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroBarriere_ImpliedVol(
		C_f,
		C_k,
		C_b,
		C_r,
		C_op,
		C_t,
		C_discount,
		C_callput,
		C_inout,
		C_updown,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroBarriere" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_discount,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	)
	{
		ADD_LOG("Local_BS_EuroBarriere_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);
	double C_discount;
	XL_readNumCell(XL_discount,C_discount,	" ARM_ERR: discount numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_inout;
	XL_GETCONVBARRIEREINOUT(XL_inout,C_inout," ARM_ERR: IN or OUT flag string expected",C_result);
	double C_updown;
	XL_GETCONVBARRIEREUPDOWN(XL_updown,C_updown," ARM_ERR: UP or DOWN flag string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroBarriere_Der(
		C_i,
		C_f,
		C_k,
		C_b,
		C_r,
		C_v,
		C_t,
		C_discount,
		C_callput,
		C_inout,
		C_updown,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroBarriere_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_discount,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	)
	{
		ADD_LOG("Local_BS_EuroBarriere_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);
	double C_discount;
	XL_readNumCell(XL_discount,C_discount,	" ARM_ERR: discount numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_inout;
	XL_GETCONVBARRIEREINOUT(XL_inout,C_inout," ARM_ERR: IN or OUT flag string expected",C_result);
	double C_updown;
	XL_GETCONVBARRIEREUPDOWN(XL_updown,C_updown," ARM_ERR: UP or DOWN flag string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroBarriere_Der2(
		C_i,
		C_j,
		C_f,
		C_k,
		C_b,
		C_r,
		C_v,
		C_t,
		C_discount,
		C_callput,
		C_inout,
		C_updown,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroBarriere_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		Double Barriere Option in Black and Sholes Model Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_r,
	LPXLOPER XL_b,
	LPXLOPER XL_callput
	)
	{
		ADD_LOG("Local_BS_EuroDoubleBarriere");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_bup;
	XL_readNumCell(XL_bup,C_bup,		" ARM_ERR: bup numeric expected",C_result);	
	double C_bdown;
	XL_readNumCell(XL_bdown,C_bdown,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,	" ARM_ERR: r numeric expected",C_result);		
	double C_b;
	XL_readNumCell(XL_b,C_b,	" ARM_ERR: b numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroDoubleBarriere(
		C_f,
		C_k,
		C_bup,
		C_bdown,
		C_v,
		C_t,
		C_r,
		C_b,
		C_callput,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroDoubleBarriere" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere_ImpliedVol(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_opt,
	LPXLOPER XL_t,
	LPXLOPER XL_r,
	LPXLOPER XL_b,
	LPXLOPER XL_callput
	)
	{
		ADD_LOG("Local_BS_EuroDoubleBarriere_ImpliedVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_bup;
	XL_readNumCell(XL_bup,C_bup,		" ARM_ERR: bup numeric expected",C_result);	
	double C_bdown;
	XL_readNumCell(XL_bdown,C_bdown,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_opt;
	XL_readNumCell(XL_opt,C_opt,	" ARM_ERR: opt numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_r;
	XL_readNumCell(XL_r,C_r,	" ARM_ERR: r numeric expected",C_result);		
	double C_b;
	XL_readNumCell(XL_b,C_b,	" ARM_ERR: b numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroDoubleBarriere_ImpliedVol(
		C_f,
		C_k,
		C_bup,
		C_bdown,
		C_opt,
		C_t,
		C_r,
		C_b,
		C_callput,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroDoubleBarriere_ImpliedVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_callput
	)
	{
		ADD_LOG("Local_BS_EuroDoubleBarriere_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_bup;
	XL_readNumCell(XL_bup,C_bup,		" ARM_ERR: bup numeric expected",C_result);	
	double C_bdown;
	XL_readNumCell(XL_bdown,C_bdown,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_EuroDoubleBarriere_Der(
		C_i,
		C_f,
		C_k,
		C_bup,
		C_bdown,
		C_v,
		C_t,
		C_callput,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroDoubleBarriere_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_callput
	)
	{
		ADD_LOG("Local_BS_EuroDoubleBarriere_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_bup;
	XL_readNumCell(XL_bup,C_bup,		" ARM_ERR: bup numeric expected",C_result);	
	double C_bdown;
	XL_readNumCell(XL_bdown,C_bdown,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: beta numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: rho numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_BS_EuroDoubleBarriere_Der2(
		C_i,
		C_j,
		C_f,
		C_k,
		C_bup,
		C_bdown,
		C_v,
		C_t,
		C_callput,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_EuroDoubleBarriere_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
////////////////////////////////////////////////////////////////////////////////////////////
///
///				Single partial barrier finishing at a date
///
////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_Start_SingleBarrier(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bendtime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_PartialTime_Start_SingleBarrier");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_barrier;
	XL_readNumCell(XL_barrier,C_barrier,		" ARM_ERR: bup numeric expected",C_result);	
	double C_rebate;
	XL_readNumCell(XL_rebate,C_rebate,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_bendtime;
	XL_readNumCell(XL_bendtime,C_bendtime,	" ARM_ERR: bendtime numeric expected",C_result);		
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_PartialTime_Start_SingleBarrier(
		C_f,
		C_k, 
		C_barrier, 
		C_rebate,
		C_v,
		C_bendtime,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_PartialTime_Start_SingleBarrier" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_Start_SingleBarrier_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bendtime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_PartialTime_Start_SingleBarrier_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_barrier;
	XL_readNumCell(XL_barrier,C_barrier,		" ARM_ERR: bup numeric expected",C_result);	
	double C_rebate;
	XL_readNumCell(XL_rebate,C_rebate,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_bendtime;
	XL_readNumCell(XL_bendtime,C_bendtime,	" ARM_ERR: bendtime numeric expected",C_result);		
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_PartialTime_Start_SingleBarrier_Der(
		C_i,
		C_f,
		C_k, 
		C_barrier, 
		C_rebate,
		C_v,
		C_bendtime,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_PartialTime_Start_SingleBarrier_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_Start_SingleBarrier_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bendtime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_PartialTime_Start_SingleBarrier_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
		double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_barrier;
	XL_readNumCell(XL_barrier,C_barrier,		" ARM_ERR: bup numeric expected",C_result);	
	double C_rebate;
	XL_readNumCell(XL_rebate,C_rebate,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_bendtime;
	XL_readNumCell(XL_bendtime,C_bendtime,	" ARM_ERR: bendtime numeric expected",C_result);		
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_PartialTime_Start_SingleBarrier_Der2(
		C_i,
		C_j,
		C_f,
		C_k, 
		C_barrier, 
		C_rebate,
		C_v,
		C_bendtime,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_PartialTime_Start_SingleBarrier_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



////////////////////////////////////////////////////////////////////////////////////////////
///
///				Single partial barrier starting at a date
///
////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_End_SingleBarrier(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_Starttime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_PartialTime_End_SingleBarrier");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_barrier;
	XL_readNumCell(XL_barrier,C_barrier,		" ARM_ERR: bup numeric expected",C_result);	
	double C_rebate;
	XL_readNumCell(XL_rebate,C_rebate,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_Starttime;
	XL_readNumCell(XL_Starttime,C_Starttime,	" ARM_ERR: bendtime numeric expected",C_result);		
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVPARTIALBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_PartialTime_Start_SingleBarrier(
		C_f,
		C_k, 
		C_barrier, 
		C_rebate,
		C_v,
		C_Starttime,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_PartialTime_End_SingleBarrier" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_End_SingleBarrier_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bstarttime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_PartialTime_End_SingleBarrier_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_barrier;
	XL_readNumCell(XL_barrier,C_barrier,		" ARM_ERR: bup numeric expected",C_result);	
	double C_rebate;
	XL_readNumCell(XL_rebate,C_rebate,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_bstarttime;
	XL_readNumCell(XL_bstarttime,C_bstarttime,	" ARM_ERR: bendtime numeric expected",C_result);		
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVPARTIALBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_PartialTime_End_SingleBarrier_Der(
		C_i,
		C_f,
		C_k, 
		C_barrier, 
		C_rebate,
		C_v,
		C_bstarttime,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_PartialTime_End_SingleBarrier_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_End_SingleBarrier_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bstarttime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_PartialTime_End_SingleBarrier_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
		/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
		double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: K numeric expected",C_result);	
	double C_barrier;
	XL_readNumCell(XL_barrier,C_barrier,		" ARM_ERR: bup numeric expected",C_result);	
	double C_rebate;
	XL_readNumCell(XL_rebate,C_rebate,		" ARM_ERR: bdown numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,	" ARM_ERR: v numeric expected",C_result);	
	double C_bstarttime;
	XL_readNumCell(XL_bstarttime,C_bstarttime,	" ARM_ERR: bendtime numeric expected",C_result);		
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVPARTIALBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_PartialTime_End_SingleBarrier_Der2(
		C_i,
		C_j,
		C_f,
		C_k, 
		C_barrier, 
		C_rebate,
		C_v,
		C_bstarttime,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_PartialTime_End_SingleBarrier_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///			Single Barrier on two assets 
///					
////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BS_SingleBarrier_2Asset(
	LPXLOPER XL_f1,
	LPXLOPER XL_k1,
	LPXLOPER XL_f2,
	LPXLOPER XL_k2, 
	LPXLOPER XL_v1, 
	LPXLOPER XL_v2,
	LPXLOPER XL_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_SingleBarrier_2Asset");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,		" ARM_ERR: f1 numeric expected",C_result);
	double C_k1;
	XL_readNumCell(XL_k1,C_k1,		" ARM_ERR: k1 numeric expected",C_result);	
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,		" ARM_ERR: f2 numeric expected",C_result);	
	double C_k2;
	XL_readNumCell(XL_k2,C_k2,		" ARM_ERR: k2 numeric expected",C_result);	
	double C_v1;
	XL_readNumCell(XL_v1,C_v1,	" ARM_ERR: v1 numeric expected",C_result);	
	double C_v2;
	XL_readNumCell(XL_v2,C_v2,	" ARM_ERR: v2 numeric expected",C_result);		
	double C_corr;
	XL_readNumCell(XL_corr,C_corr,	" ARM_ERR: corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_SingleBarrier_2Asset(
		C_f1,
		C_k1,
		C_f2,
		C_k2, 
		C_v1, 
		C_v2,
		C_corr,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_SingleBarrier_2Asset" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BS_SingleBarrier_2Asset_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f1,
	LPXLOPER XL_k1,
	LPXLOPER XL_f2,
	LPXLOPER XL_k2, 
	LPXLOPER XL_v1, 
	LPXLOPER XL_v2,
	LPXLOPER XL_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_SingleBarrier_2Asset_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,		" ARM_ERR: f1 numeric expected",C_result);
	double C_k1;
	XL_readNumCell(XL_k1,C_k1,		" ARM_ERR: k1 numeric expected",C_result);	
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,		" ARM_ERR: f2 numeric expected",C_result);	
	double C_k2;
	XL_readNumCell(XL_k2,C_k2,		" ARM_ERR: k2 numeric expected",C_result);	
	double C_v1;
	XL_readNumCell(XL_v1,C_v1,	" ARM_ERR: v1 numeric expected",C_result);	
	double C_v2;
	XL_readNumCell(XL_v2,C_v2,	" ARM_ERR: v2 numeric expected",C_result);		
	double C_corr;
	XL_readNumCell(XL_corr,C_corr,	" ARM_ERR: corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_SingleBarrier_2Asset_Der(
		C_i,
		C_f1,
		C_k1,
		C_f2,
		C_k2, 
		C_v1, 
		C_v2,
		C_corr,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_SingleBarrier_2Asset" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BS_SingleBarrier_2Asset_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f1,
	LPXLOPER XL_k1,
	LPXLOPER XL_f2,
	LPXLOPER XL_k2, 
	LPXLOPER XL_v1, 
	LPXLOPER XL_v2,
	LPXLOPER XL_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	)
	{
		ADD_LOG("Local_BS_SingleBarrier_2Asset_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,		" ARM_ERR: f1 numeric expected",C_result);
	double C_k1;
	XL_readNumCell(XL_k1,C_k1,		" ARM_ERR: k1 numeric expected",C_result);	
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,		" ARM_ERR: f2 numeric expected",C_result);	
	double C_k2;
	XL_readNumCell(XL_k2,C_k2,		" ARM_ERR: k2 numeric expected",C_result);	
	double C_v1;
	XL_readNumCell(XL_v1,C_v1,	" ARM_ERR: v1 numeric expected",C_result);	
	double C_v2;
	XL_readNumCell(XL_v2,C_v2,	" ARM_ERR: v2 numeric expected",C_result);		
	double C_corr;
	XL_readNumCell(XL_corr,C_corr,	" ARM_ERR: corr numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_BS_SingleBarrier_2Asset_Der2(
		C_i,
		C_j,
		C_f1,
		C_k1,
		C_f2,
		C_k2, 
		C_v1, 
		C_v2,
		C_corr,
		C_t,
		C_callput,
		C_optype,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_SingleBarrier_2Asset" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the bivariate function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Bivariate(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_rho,
	LPXLOPER XL_p,
	LPXLOPER XL_q
	)
	{
		ADD_LOG("Local_Bivariate");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",   C_result);
	double C_y;
	XL_readNumCell(XL_y,C_y,		" ARM_ERR: y numeric expected",   C_result);
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected", C_result);
	double C_p;
	double p_default = 0;
	XL_readNumCellWD(XL_p,C_p,	p_default,  " ARM_ERR: p numeric expected",   C_result);
	double C_q;
	double q_default = 0;
	XL_readNumCellWD(XL_q,C_q,	q_default,  " ARM_ERR: q numeric expected",   C_result);

	
	/// call the function
	long retCode=ARMLOCAL_Bivariate(
		C_x,
		C_y,
		C_rho,
		C_p,
		C_q,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bivariate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the gamma function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Gamma(
	LPXLOPER XL_x
	)
	{
		ADD_LOG("Local_Gamma");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Gamma(
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bivariate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the imcomplete beta function 
///					
////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ImcompleteBeta(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_z
	)
	{
		ADD_LOG("Local_ImcompleteBeta");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	double C_y;
	XL_readNumCell(XL_y,C_y,		" ARM_ERR: y numeric expected",C_result);
	double C_z;
	XL_readNumCell(XL_z,C_z,		" ARM_ERR: z numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_ImcompleteBeta(
		C_x,
		C_y,
		C_z,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bivariate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_InverseImcompleteBeta(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_z0,
	LPXLOPER XL_z
	)
	{
		ADD_LOG("Local_InverseImcompleteBeta");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	double C_y;
	XL_readNumCell(XL_y,C_y,		" ARM_ERR: y numeric expected",C_result);
	double C_z0;
	XL_readNumCell(XL_z0,C_z0,		" ARM_ERR: z0 numeric expected",C_result);
	double C_z;
	XL_readNumCell(XL_z,C_z,		" ARM_ERR: z numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_InverseImcompleteBeta(
		C_x,
		C_y,
		C_z0,
		C_z,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bivariate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric2F1(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_z0,
	LPXLOPER XL_z
	)
	{
		ADD_LOG("Local_Hypergeometric2F1");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: a numeric expected",C_result);
	double C_y;
	XL_readNumCell(XL_y,C_y,		" ARM_ERR: b numeric expected",C_result);
	double C_z0;
	XL_readNumCell(XL_z0,C_z0,		" ARM_ERR: c numeric expected",C_result);
	double C_z;
	XL_readNumCell(XL_z,C_z,		" ARM_ERR: z numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_Hypergeometric2F1(
		C_x,
		C_y,
		C_z0,
		C_z,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Hypergeometric2F1" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////




__declspec(dllexport) LPXLOPER WINAPI Local_RealPart_ComplexErf(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_n
	)
	{
		ADD_LOG("Local_RealPart_ComplexErf");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	double C_y;
	XL_readNumCell(XL_y,C_y,		" ARM_ERR: y numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_RealPart_ComplexErf(
		C_x,
		C_y,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RealPart_ComplexErf" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ImaginaryPart_ComplexErf(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_n
	)
	{
		ADD_LOG("Local_ImaginaryPart_ComplexErf");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	double C_y;
	XL_readNumCell(XL_y,C_y,		" ARM_ERR: y numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_ImaginaryPart_ComplexErf(
		C_x,
		C_y,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ImaginaryPart_ComplexErf" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_NormalCummulative_Inverse(
	LPXLOPER XL_x
	)
	{
		ADD_LOG("Local_NormalCummulative_Inverse");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_cdfNormal_Inv(
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NormalCummulative_Inverse" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_NormalCummulative(
	LPXLOPER XL_x
	)
	{
		ADD_LOG("Local_NormalCummulative");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_cdfNormal(
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NormalCummulative_Inverse" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		Generalized Heston  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_GHeston_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,		" ARM_ERR: sig numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_longtermV;
	XL_readNumCell(XL_longtermV,C_longtermV,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_theta;
	XL_readNumCell(XL_theta,C_theta,	" ARM_ERR: theta numeric expected",C_result);	
	double C_ksi;
	XL_readNumCell(XL_ksi,C_ksi,	" ARM_ERR: ksi numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	double C_muJ;
	XL_readNumCell(XL_muJ,C_muJ,	" ARM_ERR: muJ numeric expected",C_result);	
	double C_sigmaJ;
	XL_readNumCell(XL_sigmaJ,C_sigmaJ,	" ARM_ERR: sigmaJ numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_GHeston_VanillaOption(
		C_F,
		C_K,
		C_sig,
		C_t,
		C_longtermV,
		C_theta,
		C_ksi,
		C_rho,
		C_lambda,
		C_muJ,
		C_sigmaJ,
		C_callput,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GHeston_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_GHeston_VanillaOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,		" ARM_ERR: sig numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_longtermV;
	XL_readNumCell(XL_longtermV,C_longtermV,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_theta;
	XL_readNumCell(XL_theta,C_theta,	" ARM_ERR: theta numeric expected",C_result);	
	double C_ksi;
	XL_readNumCell(XL_ksi,C_ksi,	" ARM_ERR: ksi numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	double C_muJ;
	XL_readNumCell(XL_muJ,C_muJ,	" ARM_ERR: muJ numeric expected",C_result);	
	double C_sigmaJ;
	XL_readNumCell(XL_sigmaJ,C_sigmaJ,	" ARM_ERR: sigmaJ numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_GHeston_VanillaOption_Der(
		C_i,
		C_F,
		C_K,
		C_sig,
		C_t,
		C_longtermV,
		C_theta,
		C_ksi,
		C_rho,
		C_lambda,
		C_muJ,
		C_sigmaJ,
		C_callput,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GHeston_VanillaOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_GHeston_VanillaOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,		" ARM_ERR: sig numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_longtermV;
	XL_readNumCell(XL_longtermV,C_longtermV,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_theta;
	XL_readNumCell(XL_theta,C_theta,	" ARM_ERR: theta numeric expected",C_result);	
	double C_ksi;
	XL_readNumCell(XL_ksi,C_ksi,	" ARM_ERR: ksi numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,	" ARM_ERR: lambda numeric expected",C_result);	
	double C_muJ;
	XL_readNumCell(XL_muJ,C_muJ,	" ARM_ERR: muJ numeric expected",C_result);	
	double C_sigmaJ;
	XL_readNumCell(XL_sigmaJ,C_sigmaJ,	" ARM_ERR: sigmaJ numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_GHeston_VanillaOption_Der2(
		C_i,
		C_j,
		C_F,
		C_K,
		C_sig,
		C_t,
		C_longtermV,
		C_theta,
		C_ksi,
		C_rho,
		C_lambda,
		C_muJ,
		C_sigmaJ,
		C_callput,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GHeston_VanillaOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       CEV VanillaOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_VanillaOption(
		C_f,
		C_K,
		C_T,
		C_drift,
		C_sig,
		C_beta,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_VanillaOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_VanillaOption_Der(
		C_i,
		C_f,
		C_K,
		C_T,
		C_drift,
		C_sig,
		C_beta,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_VanillaOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_VanillaOption_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_T,
		C_drift,
		C_sig,
		C_beta,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       CEV DoubleBarrierOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_DoubleBarrierOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Bdown,
	LPXLOPER XL_Bup,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_DoubleBarrierOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);
	double C_Bdown;
	XL_readNumCell(XL_Bdown,C_Bdown,		" ARM_ERR: Bdown numeric expected",C_result);
	double C_Bup;
	XL_readNumCell(XL_Bup,C_Bup,		" ARM_ERR: Bup numeric expected",C_result);
	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_DoubleBarrierOption(
		C_f,
		C_K,
		C_T,
		C_Bdown,
		C_Bup,
		C_drift,
		C_sig,
		C_beta,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_DoubleBarrierOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_DoubleBarrierOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Bdown,
	LPXLOPER XL_Bup,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_DoubleBarrierOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);
	double C_Bdown;
	XL_readNumCell(XL_Bdown,C_Bdown,		" ARM_ERR: Bdown numeric expected",C_result);
	double C_Bup;
	XL_readNumCell(XL_Bup,C_Bup,		" ARM_ERR: Bup numeric expected",C_result);
	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_DoubleBarrierOption_Der(
		C_i,
		C_f,
		C_K,
		C_T,
		C_Bdown,
		C_Bup,
		C_drift,
		C_sig,
		C_beta,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_DoubleBarrierOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_DoubleBarrierOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Bdown,
	LPXLOPER XL_Bup,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_DoubleBarrierOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);
	double C_Bdown;
	XL_readNumCell(XL_Bdown,C_Bdown,		" ARM_ERR: Bdown numeric expected",C_result);
	double C_Bup;
	XL_readNumCell(XL_Bup,C_Bup,		" ARM_ERR: Bup numeric expected",C_result);
	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_DoubleBarrierOption_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_T,
		C_Bdown,
		C_Bup,
		C_drift,
		C_sig,
		C_beta,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_DoubleBarrierOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       CEV BarrierOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_BarrierOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_B,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,
	LPXLOPER XL_optype,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_BarrierOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);
	double C_B;
	XL_readNumCell(XL_B,C_B,		" ARM_ERR: B numeric expected",C_result);
	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_BarrierOption(
		C_f,
		C_K,
		C_T,
		C_B,
		C_drift,
		C_sig,
		C_beta,
		C_optype,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_BarrierOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_BarrierOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_B,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_optype,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_BarrierOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);
	double C_B;
	XL_readNumCell(XL_B,C_B,		" ARM_ERR: B numeric expected",C_result);
	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_BarrierOption_Der(
		C_i,
		C_f,
		C_K,
		C_T,
		C_B,
		C_drift,
		C_sig,
		C_beta,
		C_optype,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_BarrierOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_BarrierOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_B,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,
	LPXLOPER XL_optype,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_CEV_BarrierOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);
	double C_B;
	XL_readNumCell(XL_B,C_B,		" ARM_ERR: B numeric expected",C_result);
	double C_drift;
	XL_readNumCell(XL_drift,C_drift, " ARM_ERR: drift numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);
	double C_optype;
	XL_GETCONVBARRIERETYPE(XL_optype,C_optype," ARM_ERR: optype string expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_CEV_BarrierOption_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_T,
		C_B,
		C_drift,
		C_sig,
		C_beta,
		C_optype,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_BarrierOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadDigitalOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a1;
	XL_readNumCell(XL_a1,C_a1,	" ARM_ERR: a1 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);

	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadDigitalOption(
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a1,
		C_a2,
		C_a3,
		C_t,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadDigitalOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadDigitalOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a1;
	XL_readNumCell(XL_a1,C_a1,	" ARM_ERR: a1 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);

	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadDigitalOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a1,
		C_a2,
		C_a3,
		C_t,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadDigitalOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadDigitalOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a1;
	XL_readNumCell(XL_a1,C_a1,	" ARM_ERR: a1 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);

	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadDigitalOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a1,
		C_a2,
		C_a3,
		C_t,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadDigitalOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option2 Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption2(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_b0,	
	LPXLOPER XL_b1,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadDigitalOption2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);
	double C_b0;
	XL_readNumCell(XL_b0,C_b0,	" ARM_ERR: b0 numeric expected",C_result);
	double C_b1;
	XL_readNumCell(XL_b1,C_b1,	" ARM_ERR: b1 numeric expected",C_result);
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	

	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadDigitalOption2(
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a2,
		C_a3,
		C_b0,
		C_b1,
		C_t,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadDigitalOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption2_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_b0,	
	LPXLOPER XL_b1,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadDigitalOption2_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);
	double C_b0;
	XL_readNumCell(XL_b0,C_b0,	" ARM_ERR: b0 numeric expected",C_result);
	double C_b1;
	XL_readNumCell(XL_b1,C_b1,	" ARM_ERR: b1 numeric expected",C_result);
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadDigitalOption2_Der(
		C_i,
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a2,
		C_a3,
		C_b0,
		C_b1,
		C_t,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadDigitalOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption2_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_b0,	
	LPXLOPER XL_b1,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadDigitalOption2_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	
	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);
	double C_b0;
	XL_readNumCell(XL_b0,C_b0,	" ARM_ERR: b0 numeric expected",C_result);
	double C_b1;
	XL_readNumCell(XL_b1,C_b1,	" ARM_ERR: b1 numeric expected",C_result);
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadDigitalOption2_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a2,
		C_a3,
		C_b0,
		C_b1,
		C_t,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadDigitalOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a1;
	XL_readNumCell(XL_a1,C_a1,	" ARM_ERR: a1 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);

	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadOption(
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a1,
		C_a2,
		C_a3,
		C_t,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a1;
	XL_readNumCell(XL_a1,C_a1,	" ARM_ERR: a1 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);

	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadOption_Der(
		C_i,
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a1,
		C_a2,
		C_a3,
		C_t,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_LogNormal_TriSpreadOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);	
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S1 numeric expected",C_result);
	double C_S3;
	XL_readNumCell(XL_S3,C_S3,		" ARM_ERR: S1 numeric expected",C_result);

	double C_sig1;
	XL_readNumCell(XL_sig1,C_sig1,	" ARM_ERR: sig1 numeric expected",C_result);	
	double C_sig2;
	XL_readNumCell(XL_sig2,C_sig2,	" ARM_ERR: sig2 numeric expected",C_result);	
	double C_sig3;
	XL_readNumCell(XL_sig3,C_sig3,	" ARM_ERR: sig3 numeric expected",C_result);
	
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,	" ARM_ERR: rho12 numeric expected",C_result);	
		double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,	" ARM_ERR: rho numeric expected",C_result);	
		double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,	" ARM_ERR: rho numeric expected",C_result);	

	double C_mu1;
	XL_readNumCell(XL_mu1,C_mu1,	" ARM_ERR: mu1 numeric expected",C_result);	
	double C_mu2;
	XL_readNumCell(XL_mu2,C_mu2,	" ARM_ERR: mu2 numeric expected",C_result);	
	double C_mu3;
	XL_readNumCell(XL_mu3,C_mu3,	" ARM_ERR: mu3 numeric expected",C_result);	

	double C_a0;
	XL_readNumCell(XL_a0,C_a0,	" ARM_ERR: a0 numeric expected",C_result);
		double C_a1;
	XL_readNumCell(XL_a1,C_a1,	" ARM_ERR: a1 numeric expected",C_result);
		double C_a2;
	XL_readNumCell(XL_a2,C_a2,	" ARM_ERR: a2 numeric expected",C_result);
		double C_a3;
	XL_readNumCell(XL_a3,C_a3,	" ARM_ERR: a3 numeric expected",C_result);

	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: t numeric expected",C_result);	

	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_LogNormal_TriSpreadOption_Der2(
		C_i,
		C_j,
		C_S1,
		C_S2,
		C_S3,
		C_sig1,
		C_sig2,
		C_sig3,
		C_rho12,
		C_rho13,
		C_rho23,
		C_mu1,
		C_mu2,
		C_mu3,
		C_a0,
		C_a1,
		C_a2,
		C_a3,
		C_t,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LogNormal_TriSpreadOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Asian Lognormal vanilla option (geman yor formula)
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Lognormal_Asian_VanillaOption(
	LPXLOPER XL_S,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_alpha,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_Lognormal_Asian_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S;
	XL_readNumCell(XL_S,C_S,		" ARM_ERR: S numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,	" ARM_ERR: theta numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Lognormal_Asian_VanillaOption(
		C_S,
		C_k,
		C_T,
		C_r,
		C_v,
		C_alpha,
		C_callput,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Lognormal_Asian_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Lognormal_Asian_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_alpha,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_Lognormal_Asian_VanillaOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_S;
	XL_readNumCell(XL_S,C_S,		" ARM_ERR: S numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,	" ARM_ERR: theta numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Lognormal_Asian_VanillaOption_Der(
		C_i,
		C_S,
		C_k,
		C_T,
		C_r,
		C_v,
		C_alpha,
		C_callput,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Lognormal_Asian_VanillaOption_Der" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Lognormal_Asian_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_alpha,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_Lognormal_Asian_VanillaOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: i numeric expected",C_result);
	double C_S;
	XL_readNumCell(XL_S,C_S,		" ARM_ERR: S numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,		" ARM_ERR: k numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_r;
	XL_readNumCell(XL_r,C_r,		" ARM_ERR: r numeric expected",C_result);	
	double C_v;
	XL_readNumCell(XL_v,C_v,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,	" ARM_ERR: theta numeric expected",C_result);		
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Lognormal_Asian_VanillaOption_Der2(
		C_i,
		C_j,
		C_S,
		C_k,
		C_T,
		C_r,
		C_v,
		C_alpha,
		C_callput,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Lognormal_Asian_VanillaOption_Der2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GaussianIntegrals_Legendre_Coeffs(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_n)
{
	ADD_LOG("Local_GaussianIntegrals_Legendre_Coeffs");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_a;
		XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);
		double C_b;
		XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);

		double C_n;
		XL_readNumCell(XL_n,C_n,		" ARM_ERR: n numeric expected",C_result);

		
		/// call the function
		long retCode=ARMLOCAL_GaussianIntegrals_Legendre_Coeffs(
			C_a,
			C_b,
			C_n,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = (int)C_n;
			int nbCols = 2;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GaussianIntegrals_Legendre_Coeffs" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GaussianIntegrals_Hermite_Coeffs(
	LPXLOPER XL_n )
{
	ADD_LOG("Local_GaussianIntegrals_Hermite_Coeffs");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_n;
		XL_readNumCell(XL_n,C_n,		" ARM_ERR: n numeric expected",C_result);
		
		/// call the function
		long retCode=ARMLOCAL_GaussianIntegrals_Hermite_Coeffs(
			C_n,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = (int)C_n;
			int nbCols = 2;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GaussianIntegrals_Hermite_Coeffs" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GaussianIntegrals_Laguerre_Coeffs(
	LPXLOPER XL_a,
	LPXLOPER XL_n)
{
	ADD_LOG("Local_GaussianIntegrals_Laguerre_Coeffs");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_n;
		XL_readNumCell(XL_n,C_n,		" ARM_ERR: n numeric expected",C_result);
		double C_a;
		XL_readNumCell(XL_a,C_a,		" ARM_ERR: n numeric expected",C_result);
		
		/// call the function
		long retCode=ARMLOCAL_GaussianIntegrals_Laguerre_Coeffs(
			C_a,
			C_n,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = (int)C_n;
			int nbCols = 2;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GaussianIntegrals_Laguerre_Coeffs" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// function to calibrate in an FX option the strike from the modified delta
/////////////////////////////////////////////////////////////
LPXLOPER Local_ConvertFXOptionToStrike(
    LPXLOPER XL_TargetDelta,
	LPXLOPER XL_Fwd,
	LPXLOPER XL_TotalVol,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_InitValue )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		double C_TargetPrice;
		XL_readNumCell(XL_TargetDelta,C_TargetPrice,	" ARM_ERR: Target Price numeric expected",C_result);
		double C_Fwd;
		XL_readNumCell(XL_Fwd,C_Fwd,					" ARM_ERR: Fwd numeric expected",C_result);
		double C_TotalVol;
		XL_readNumCell(XL_TotalVol,C_TotalVol,			" ARM_ERR: TotalVol numeric expected",C_result);
		double C_callput;
		XL_GETCONVCALLORPUT(XL_CallPut,C_callput,		" ARM_ERR: call or put string expected",C_result);
		double C_InitValue;
		XL_readNumCell(XL_InitValue,C_InitValue,		" ARM_ERR: Init Value numeric expected",C_result);

		/// a function with a context
		long retCode=ARMLOCAL_ConvertFXOptionToStrike( 
			C_TargetPrice,
			C_Fwd,
			C_TotalVol,
			C_callput,
			C_InitValue,
			C_result );

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ConvertFXOptionToStrike" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/////////////////////////////////////////////////////////////////////////////////////
///
///		Generalized Heston  Vanilla option Formula With Vector of Model
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption_ModelVector(
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_Interpolation_Method,
	LPXLOPER XL_Maturity_Vec,
	LPXLOPER XL_F_Vec,
	LPXLOPER XL_InitialVol_Vec,
	LPXLOPER XL_longtermV_Vec,
	LPXLOPER XL_theta_Vec,
	LPXLOPER XL_ksi_Vec,
	LPXLOPER XL_rho_Vec,
	LPXLOPER XL_lambda_Vec,
	LPXLOPER XL_muJ_Vec,
	LPXLOPER XL_sigmaJ_Vec,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_GHeston_VanillaOption_ModelVector");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_Interpolation_Method;
	XL_GETCONVHESTONINTERPOMETH(XL_Interpolation_Method,C_Interpolation_Method," ARM_ERR: Interpolation Method expected",C_result);

	vector<double> C_Maturity_Vec;
	XL_readNumVector(XL_Maturity_Vec,C_Maturity_Vec," ARM_ERR: Maturity Vector: array of numeric expected",C_result);

	vector<double> C_F_Vec;
	XL_readNumVector(XL_F_Vec,C_F_Vec,		" ARM_ERR: F Vector: array of  numeric expected",C_result);	

	vector<double> C_InitialVol_Vec;
	XL_readNumVector(XL_InitialVol_Vec,C_InitialVol_Vec,		" ARM_ERR: InitialVol Vector: array of  numeric expected",C_result);	

	vector<double> C_longtermV_Vec;
	XL_readNumVector(XL_longtermV_Vec,C_longtermV_Vec,		" ARM_ERR: longtermV Vector: array of  numeric expected",C_result);	

	vector<double> C_theta_Vec;
	XL_readNumVector(XL_theta_Vec,C_theta_Vec,	" ARM_ERR: theta Vector: array of  numeric expected",C_result);	

	vector<double> C_ksi_Vec;
	XL_readNumVector(XL_ksi_Vec,C_ksi_Vec,	" ARM_ERR: ksi Vector: array of  numeric expected",C_result);	

	vector<double> C_rho_Vec;
	XL_readNumVector(XL_rho_Vec,C_rho_Vec,	" ARM_ERR: rho Vector: array of  numeric expected",C_result);	

	vector<double> C_lambda_Vec;
	XL_readNumVector(XL_lambda_Vec,C_lambda_Vec,	" ARM_ERR: lambda Vector: array of  numeric expected",C_result);	

	vector<double> C_muJ_Vec;
	XL_readNumVector(XL_muJ_Vec,C_muJ_Vec,	" ARM_ERR: muJ Vector: array of  numeric expected",C_result);	

	vector<double> C_sigmaJ_Vec;
	XL_readNumVector(XL_sigmaJ_Vec,C_sigmaJ_Vec,	" ARM_ERR: sigmaJ Vector: array of  numeric expected",C_result);	

	double C_nb;
	double n_default=80;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_GHeston_VanillaOption_ModelVector(
		C_K,
		C_t,
		C_callput,
		C_Interpolation_Method,
		C_Maturity_Vec,
		C_F_Vec,
		C_InitialVol_Vec,
		C_longtermV_Vec,
		C_theta_Vec,
		C_ksi_Vec,
		C_rho_Vec,
		C_lambda_Vec,
		C_muJ_Vec,
		C_sigmaJ_Vec,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GHeston_VanillaOption_ModelVector" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_Implicit_Volatility_ModelVector(
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_Interpolation_Method,
	LPXLOPER XL_Maturity_Vec,
	LPXLOPER XL_F_Vec,
	LPXLOPER XL_InitialVol_Vec,
	LPXLOPER XL_longtermV_Vec,
	LPXLOPER XL_theta_Vec,
	LPXLOPER XL_ksi_Vec,
	LPXLOPER XL_rho_Vec,
	LPXLOPER XL_lambda_Vec,
	LPXLOPER XL_muJ_Vec,
	LPXLOPER XL_sigmaJ_Vec,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_GHeston_Implicit_Volatility_ModelVector");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_Interpolation_Method;
	XL_GETCONVHESTONINTERPOMETH(XL_Interpolation_Method,C_Interpolation_Method," ARM_ERR: Interpolation Method expected",C_result);

	vector<double> C_Maturity_Vec;
	XL_readNumVector(XL_Maturity_Vec,C_Maturity_Vec," ARM_ERR: Maturity Vector: array of numeric expected",C_result);

	vector<double> C_F_Vec;
	XL_readNumVector(XL_F_Vec,C_F_Vec,		" ARM_ERR: F Vector: array of  numeric expected",C_result);	

	vector<double> C_InitialVol_Vec;
	XL_readNumVector(XL_InitialVol_Vec,C_InitialVol_Vec,		" ARM_ERR: InitialVol Vector: array of  numeric expected",C_result);	

	vector<double> C_longtermV_Vec;
	XL_readNumVector(XL_longtermV_Vec,C_longtermV_Vec,		" ARM_ERR: longtermV Vector: array of  numeric expected",C_result);	

	vector<double> C_theta_Vec;
	XL_readNumVector(XL_theta_Vec,C_theta_Vec,	" ARM_ERR: theta Vector: array of  numeric expected",C_result);	

	vector<double> C_ksi_Vec;
	XL_readNumVector(XL_ksi_Vec,C_ksi_Vec,	" ARM_ERR: ksi Vector: array of  numeric expected",C_result);	

	vector<double> C_rho_Vec;
	XL_readNumVector(XL_rho_Vec,C_rho_Vec,	" ARM_ERR: rho Vector: array of  numeric expected",C_result);	

	vector<double> C_lambda_Vec;
	XL_readNumVector(XL_lambda_Vec,C_lambda_Vec,	" ARM_ERR: lambda Vector: array of  numeric expected",C_result);	

	vector<double> C_muJ_Vec;
	XL_readNumVector(XL_muJ_Vec,C_muJ_Vec,	" ARM_ERR: muJ Vector: array of  numeric expected",C_result);	

	vector<double> C_sigmaJ_Vec;
	XL_readNumVector(XL_sigmaJ_Vec,C_sigmaJ_Vec,	" ARM_ERR: sigmaJ Vector: array of  numeric expected",C_result);	

	double C_nb;
	double n_default=80;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_GHeston_Implicit_Volatility_ModelVector(
		C_K,
		C_t,
		C_Interpolation_Method,
		C_Maturity_Vec,
		C_F_Vec,
		C_InitialVol_Vec,
		C_longtermV_Vec,
		C_theta_Vec,
		C_ksi_Vec,
		C_rho_Vec,
		C_lambda_Vec,
		C_muJ_Vec,
		C_sigmaJ_Vec,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GHeston_VanillaOption_ModelVector" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		Calibration of SABR
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate_BetaFixed(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_beta,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_Weight_Vec,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh

	)
{
	ADD_LOG("Local_SABR_Model_Calibrate_BetaFixed");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
		double C_f;
		XL_readNumCell(XL_f,C_f,		" ARM_ERR: f Value numeric expected",C_result);
		double C_beta;
		XL_readNumCell(XL_beta,C_beta,		" ARM_ERR: beta Value numeric expected",C_result);
		double C_t;
		XL_readNumCell(XL_t,C_t,		" ARM_ERR: t Value numeric expected",C_result);

		double C_flag;
		XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
		double C_nb;
		XL_readNumCell(XL_nbsteps, C_nb,		 " ARM_ERR: nb: numeric expected",C_result);
	
		double C_algorithm;
		XL_GETCONVOPTIMIZATIONFLAG(XL_algorithm,C_algorithm," ARM_ERR: Algorithm flag string expected",C_result);
		double C_alpha0;
		XL_readNumCell(XL_alpha0, C_alpha0,		 " ARM_ERR: alpha0: numeric expected",C_result);
		double C_rho0;
		XL_readNumCell(XL_rho0, C_rho0,		 " ARM_ERR: beta0: numeric expected",C_result);
		double C_nu0;
		XL_readNumCell(XL_nu0, C_nu0,		 " ARM_ERR: nu0: numeric expected",C_result);
		vector<double> C_Weight_Vec;
	
		int sizevec=C_ImpVol_Vec.size();
		vector<double> C_defaultVector(sizevec,1.0);
		XL_readNumVectorWD(XL_Weight_Vec,C_Weight_Vec,C_defaultVector," ARM_ERR: Weight Vector: array of numeric expected",C_result);
		double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);
		/// call the function
		long retCode=ARMLOCAL_SABR_Model_Calibrate_FixedBeta(
			C_K_Vec,
			C_ImpVol_Vec,
			C_Weight_Vec,
			C_f,
			C_beta,
			C_t,
			C_flag,
			C_Alpha_Exp,
			C_Alpha_Tanh,
			C_Kb_Tanh,
			C_nb,
			C_algorithm,
			C_alpha0,
			C_rho0,
			C_nu0,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate_BetaFixed_Linked(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_beta,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_alphap,
	LPXLOPER XL_rhop,
	LPXLOPER XL_nup,
	LPXLOPER XL_rweight_alpha,
	LPXLOPER XL_rweight_rho,
	LPXLOPER XL_rweight_nu,
	LPXLOPER XL_Weight_Vec
	)
{
	ADD_LOG("Local_SABR_Model_Calibrate_BetaFixed_Linked");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
		double C_f;
		XL_readNumCell(XL_f,C_f,		" ARM_ERR: f Value numeric expected",C_result);
		double C_beta;
		XL_readNumCell(XL_beta,C_beta,		" ARM_ERR: beta Value numeric expected",C_result);
		double C_t;
		XL_readNumCell(XL_t,C_t,		" ARM_ERR: t Value numeric expected",C_result);

		double C_flag;
		XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
		double C_nb;
		XL_readNumCell(XL_nbsteps, C_nb,		 " ARM_ERR: nb: numeric expected",C_result);
	
		double C_algorithm;
		XL_GETCONVOPTIMIZATIONFLAG(XL_algorithm,C_algorithm," ARM_ERR: Algorithm flag string expected",C_result);
		double C_alpha0;
		XL_readNumCell(XL_alpha0, C_alpha0,		 " ARM_ERR: alpha0: numeric expected",C_result);
		double C_rho0;
		XL_readNumCell(XL_rho0, C_rho0,		 " ARM_ERR: beta0: numeric expected",C_result);
		double C_nu0;
		XL_readNumCell(XL_nu0, C_nu0,		 " ARM_ERR: nu0: numeric expected",C_result);
		double C_alphap;
		XL_readNumCell(XL_alphap, C_alphap,		 " ARM_ERR: alphap: numeric expected",C_result);
		double C_rhop;
		XL_readNumCell(XL_rhop, C_rhop,		 " ARM_ERR: betap: numeric expected",C_result);
		double C_nup;
		XL_readNumCell(XL_nup, C_nup,		 " ARM_ERR: nup: numeric expected",C_result);
		double C_rweight_alpha;
		XL_readNumCell(XL_rweight_alpha, C_rweight_alpha,		 " ARM_ERR: rweight_alpha: numeric expected",C_result);
		double C_rweight_rho;
		XL_readNumCell(XL_rweight_rho, C_rweight_rho,		 " ARM_ERR: rweight_rho: numeric expected",C_result);
		double C_rweight_nu;
		XL_readNumCell(XL_rweight_nu, C_rweight_nu,		 " ARM_ERR: rweight_nu: numeric expected",C_result);
		vector<double> C_Weight_Vec;
		XL_readNumVector(XL_Weight_Vec,C_Weight_Vec," ARM_ERR: Weight Vector: array of numeric expected",C_result);
	
	
		/// call the function
		long retCode=ARMLOCAL_SABR_Model_Calibrate_FixedBeta_Linked(
			C_K_Vec,
			C_ImpVol_Vec,
			C_Weight_Vec,
			C_f,
			C_beta,
			C_t,
			C_flag,
			C_nb,
			C_algorithm,
			C_alpha0,
			C_rho0,
			C_nu0,
			C_alphap,
			C_rhop,
			C_nup,
			C_rweight_alpha,
			C_rweight_rho,
			C_rweight_nu,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_beta0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_Weight_Vec
	)
{
	ADD_LOG("Local_SABR_Model_Calibrate");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
			double C_f;
		XL_readNumCell(XL_f,C_f,		" ARM_ERR: f Value numeric expected",C_result);

		double C_t;
		XL_readNumCell(XL_t,C_t,		" ARM_ERR: t Value numeric expected",C_result);

		double C_flag;
		XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
		double C_nb;
		XL_readNumCell(XL_nbsteps, C_nb,		 " ARM_ERR: nb: numeric expected",C_result);
	
		double C_algorithm;
		XL_GETCONVOPTIMIZATIONFLAG(XL_algorithm,C_algorithm," ARM_ERR: Algorithm flag string expected",C_result);
		double C_alpha0;
		XL_readNumCell(XL_alpha0, C_alpha0,		 " ARM_ERR: alpha0: numeric expected",C_result);
		double C_beta0;
		XL_readNumCell(XL_beta0, C_beta0,		 " ARM_ERR: beta0: numeric expected",C_result);
		double C_rho0;
		XL_readNumCell(XL_rho0, C_rho0,		 " ARM_ERR: rho0: numeric expected",C_result);
		double C_nu0;
		XL_readNumCell(XL_nu0, C_nu0,		 " ARM_ERR: nu0: numeric expected",C_result);
		vector<double> C_Weight_Vec;
		XL_readNumVector(XL_Weight_Vec,C_Weight_Vec," ARM_ERR: Weight Vector: array of numeric expected",C_result);
	
	
		/// call the function
		long retCode=ARMLOCAL_SABR_Model_Calibrate(
			C_K_Vec,
			C_ImpVol_Vec,
			C_Weight_Vec,
			C_f,
			C_t,
			C_flag,
			C_nb,
			C_algorithm,
			C_alpha0,
			C_beta0,
			C_rho0,
			C_nu0,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate_WithForward(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_beta0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_f0,
	LPXLOPER XL_Weight_Vec
	)
{
	ADD_LOG("Local_SABR_Model_Calibrate_WithForward");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
		double C_t;
		XL_readNumCell(XL_t,C_t,		" ARM_ERR: t Value numeric expected",C_result);

		double C_flag;
		XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
		double C_nb;
		XL_readNumCell(XL_nbsteps, C_nb,		 " ARM_ERR: nb: numeric expected",C_result);
		double C_algorithm;
		XL_GETCONVOPTIMIZATIONFLAG(XL_algorithm,C_algorithm," ARM_ERR: Algorithm flag string expected",C_result);
		double C_alpha0;
		XL_readNumCell(XL_alpha0, C_alpha0,		 " ARM_ERR: alpha0: numeric expected",C_result);
		double C_beta0;
		XL_readNumCell(XL_beta0, C_beta0,		 " ARM_ERR: beta0: numeric expected",C_result);
		double C_rho0;
		XL_readNumCell(XL_rho0, C_rho0,		 " ARM_ERR: rho0: numeric expected",C_result);
		double C_nu0;
		XL_readNumCell(XL_nu0, C_nu0,		 " ARM_ERR: nu0: numeric expected",C_result);
		double C_f0;
		XL_readNumCell(XL_f0, C_f0,		 " ARM_ERR: f0: numeric expected",C_result);
		vector<double> C_Weight_Vec;
		XL_readNumVector(XL_Weight_Vec,C_Weight_Vec," ARM_ERR: Weight Vector: array of numeric expected",C_result);
	
	
		
		/// call the function
		long retCode=ARMLOCAL_SABR_Model_Calibrate_WForward(
			C_K_Vec,
			C_ImpVol_Vec,
			C_Weight_Vec,
			C_t,
			C_flag,
			C_nb,
			C_algorithm,
			C_alpha0,
			C_beta0,
			C_rho0,
			C_nu0,
			C_f0,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Calibrate_BetaFixedToOne(
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_fwd,
	LPXLOPER XL_mat,
	LPXLOPER XL_atmvol,
	LPXLOPER XL_Weight_Vec
	)
{
	ADD_LOG("Local_SABR_Calibrate_BetaFixedToOne");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
		double C_t;
		XL_readNumCell(XL_mat,C_t,		" ARM_ERR: t Value numeric expected",C_result);

		double C_fwd;
		XL_readNumCell(XL_fwd,C_fwd,		" ARM_ERR: fwd Value numeric expected",C_result);

		double C_atmvol,defaultatm = 0.;
		XL_readNumCellWD(XL_atmvol,C_atmvol,defaultatm, "ARM_ERR : atm vol : numeric expected",C_result);

		vector<double> C_Weight_Vec, DefaultVec;
		XL_readNumVectorWD(XL_Weight_Vec,C_Weight_Vec,DefaultVec," ARM_ERR: Weight Vector: array of numeric expected",C_result);
	
	
		/// call the function
		long retCode=ARMLOCAL_SABR_Calibrate_BetaFixedToOne(
			C_K_Vec,
			C_ImpVol_Vec,
			C_Weight_Vec,
			C_fwd,
			C_t,
			C_atmvol,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/////////////////////////////////////////////////////////////////////////////////////
///
///		Calibration of GHeston
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_Model_Calibrate_Total(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_t,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_V0	,
	LPXLOPER XL_omega,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi	,
	LPXLOPER XL_rho	,
	LPXLOPER XL_muJ	,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_lambda
	)
{
	ADD_LOG("Local_GHeston_Model_Calibrate_Total");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
		double C_f;
		XL_readNumCell(XL_f, C_f,		 " ARM_ERR: f: numeric expected",C_result);
		double C_t;
		XL_readNumCell(XL_t,C_t,		" ARM_ERR: t Value numeric expected",C_result);
		double C_nb;
		XL_readNumCell(XL_nbsteps, C_nb,		 " ARM_ERR: nb: numeric expected",C_result);
		double C_algorithm;
		XL_GETCONVOPTIMIZATIONFLAG(XL_algorithm,C_algorithm," ARM_ERR: Algorithm flag string expected",C_result);
		double C_V0;
		XL_readNumCell(XL_V0, C_V0,		 " ARM_ERR: V0: numeric expected",C_result);
		double C_omega;
		XL_readNumCell(XL_omega, C_omega,		 " ARM_ERR: omega: numeric expected",C_result);
		double C_theta;
		XL_readNumCell(XL_theta, C_theta,		 " ARM_ERR: theta: numeric expected",C_result);
		double C_ksi;
		XL_readNumCell(XL_ksi, C_ksi,		 " ARM_ERR: ksi: numeric expected",C_result);
			double C_rho;
		XL_readNumCell(XL_rho, C_rho,		 " ARM_ERR: rho: numeric expected",C_result);
			double C_muJ;
		XL_readNumCell(XL_muJ, C_muJ,		 " ARM_ERR: muJ: numeric expected",C_result);
			double C_sigmaJ;
		XL_readNumCell(XL_sigmaJ, C_sigmaJ,		 " ARM_ERR: sigmaJ: numeric expected",C_result);
			double C_lambda;
		XL_readNumCell(XL_lambda, C_lambda,		 " ARM_ERR: Lambda: numeric expected",C_result);
	
	
		
		/// call the function
		long retCode=ARMLOCAL_GHeston_Model_Calibrate_Total(
			C_K_Vec,
			C_ImpVol_Vec,
			C_f,
			C_t,
			C_nb,
			C_algorithm,
			C_V0,
			C_omega,
			C_theta,
			C_ksi,
			C_rho,
			C_muJ,
			C_sigmaJ,
			C_lambda,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 9;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GHeston_Model_Calibrate_Total" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_Model_Calibrate_NoJump(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_t,
	LPXLOPER XL_muJ0	,
	LPXLOPER XL_sigmaJ0,
	LPXLOPER XL_lambda0,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_V0	,
	LPXLOPER XL_omega,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi	,
	LPXLOPER XL_rho	
	)
{
	ADD_LOG("Local_GHeston_Model_Calibrate_NoJump");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_K_Vec;
		XL_readNumVector(XL_K_Vec,C_K_Vec," ARM_ERR: K Vector: array of numeric expected",C_result);
		
		vector<double> C_ImpVol_Vec;
		XL_readNumVector(XL_ImpVol_Vec,C_ImpVol_Vec," ARM_ERR: ImpVol Vector: array of numeric expected",C_result);
		double C_f;
		XL_readNumCell(XL_f, C_f,		 " ARM_ERR: f: numeric expected",C_result);
		double C_t;
		XL_readNumCell(XL_t,C_t,		" ARM_ERR: t Value numeric expected",C_result);
		double C_nb;
		XL_readNumCell(XL_nbsteps, C_nb,		 " ARM_ERR: nb: numeric expected",C_result);
		double C_algorithm;
		XL_GETCONVOPTIMIZATIONFLAG(XL_algorithm,C_algorithm," ARM_ERR: Algorithm flag string expected",C_result);
		double C_V0;
		XL_readNumCell(XL_V0, C_V0,		 " ARM_ERR: V0: numeric expected",C_result);
		double C_omega;
		XL_readNumCell(XL_omega, C_omega,		 " ARM_ERR: omega: numeric expected",C_result);
		double C_theta;
		XL_readNumCell(XL_theta, C_theta,		 " ARM_ERR: theta: numeric expected",C_result);
		double C_ksi;
		XL_readNumCell(XL_ksi, C_ksi,		 " ARM_ERR: ksi: numeric expected",C_result);
			double C_rho;
		XL_readNumCell(XL_rho, C_rho,		 " ARM_ERR: rho: numeric expected",C_result);
			double C_muJ0;
		XL_readNumCell(XL_muJ0, C_muJ0,		 " ARM_ERR: muJ0: numeric expected",C_result);
			double C_sigmaJ0;
		XL_readNumCell(XL_sigmaJ0, C_sigmaJ0,		 " ARM_ERR: sigmaJ0: numeric expected",C_result);
			double C_lambda0;
		XL_readNumCell(XL_lambda0, C_lambda0,		 " ARM_ERR: Lambda: numeric expected",C_result);
	
	
		
		/// call the function
		long retCode=ARMLOCAL_GHeston_Model_Calibrate_NoJump(
			C_K_Vec,
			C_ImpVol_Vec,
			C_f,
			C_t,
			C_nb,
			C_algorithm,
			C_V0,
			C_omega,
			C_theta,
			C_ksi,
			C_rho,
			C_muJ0,
			C_sigmaJ0,
			C_lambda0,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 8;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR + Student law
///
/////////////////////////////////////////////////////////////////////////////////////




__declspec(dllexport) LPXLOPER WINAPI Local_Student_SABR_Power_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_Parameters_Corr,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec
	)
{
	ADD_LOG("Local_Student_SABR_Power_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	vector<double> C_Parameters_Corr;
	XL_readNumVector(XL_Parameters_Corr,C_Parameters_Corr," ARM_ERR: Parameters Corr Vector : array of 2 numeric expected",C_result);
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_a10;
	XL_readNumCell(XL_a10,C_a10,		" ARM_ERR: a10 numeric expected",C_result);
	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);	
	vector<double> C_Parameters_vec;
	XL_readNumVector(XL_Parameters_vec,C_Parameters_vec," ARM_ERR: Parameters Vector : array of numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_Student_SABR_Power_SpreadOption(
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_Parameters_Corr,
		C_t,
		C_flag,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_Parameters_vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Gaussian_SABR_Power_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_Student_SABR_Power_Digital_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh
	)
{
	ADD_LOG("Local_Student_SABR_Power_Digital_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,		" ARM_ERR: S1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,		" ARM_ERR: S2 numeric expected",C_result);	
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,	" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,	" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,	" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,	" ARM_ERR: nu1 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,	" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,	" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,	" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,	" ARM_ERR: nu2 numeric expected",C_result);	
	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_copula_degre;
	XL_readNumCell(XL_copula_degre,C_copula_degre,		" ARM_ERR: copula_degre numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: t numeric expected",C_result);
	double C_flag;
	XL_GETCONVSABRFLAG(XL_flag,C_flag," ARM_ERR: SABR flag string expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default=120;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);
	double C_alpha_exp;
	double alpha_exp_default=0.01;
    XL_readNumCellWD(XL_alpha_exp, C_alpha_exp, alpha_exp_default, " ARM_ERR: n: numeric expected",C_result);
	double C_alpha_tanh;
	double alpha_tanh_default=1.5;
    XL_readNumCellWD(XL_alpha_tanh, C_alpha_tanh, alpha_tanh_default, " ARM_ERR: n: numeric expected",C_result);
	double C_kb_tanh;
	double kb_tanh_default=0.02;
    XL_readNumCellWD(XL_kb_tanh, C_kb_tanh, kb_tanh_default, " ARM_ERR: kb_tanh: numeric expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_Student_SABR_Power_Digital_SpreadOption(
		C_S1,
		C_S2,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_copula_corr,
		C_copula_degre,
		C_t,
		C_flag,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_n,C_alpha_exp,C_alpha_tanh, C_kb_tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_SABR_Power_Digital_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//////////////////////////////////////////////////////////////////////////////////


/// Mepi Vanilla Option 


/////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_STOBS(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nz,
	LPXLOPER XL_nh)
{
	ADD_LOG("Local_Mepi_VanillaOption_STOBS");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs	
		double C_date;
	    XL_readNumCell(XL_date,C_date,		" ARM_ERR: date numeric expected",C_result);	
	    double C_P0;
	    XL_readNumCell(XL_P0,C_P0,		" ARM_ERR: P0 numeric expected",C_result);	
	    double C_K;
	    XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	    double C_T;
	    XL_readNumCell(XL_T,C_T,	            " ARM_ERR: T numeric expected",C_result);	
		

		long C_Zero_Curve_id;		
		XL_GETOBJID( XL_ZeroCurve, C_Zero_Curve_id,	" ARM_ERR: Zero Curve: Object expected",		C_result);
	
		double C_b;
		XL_readNumCell(XL_b,C_b,	        " ARM_ERR: b numeric expected",C_result);
		double C_YearlyFees;
		XL_readNumCell(XL_YearlyFees,C_YearlyFees,	        " ARM_ERR: YearlyFees numeric expected",C_result);
		
		
		long C_SABR_Model_id;		
		XL_GETOBJID( XL_SABR, C_SABR_Model_id,	" ARM_ERR: SABR Model: Object expected",		C_result);
		
	
		double C_minExp;
		XL_readNumCell(XL_minExp,C_minExp,	        " ARM_ERR: minExp numeric expected",C_result);
		double C_maxExp;
		XL_readNumCell(XL_maxExp,C_maxExp,	        " ARM_ERR: maxExp numeric expected",C_result);
		double C_riskFac;
		XL_readNumCell(XL_riskFac,C_riskFac,	        " ARM_ERR: riskFac numeric expected",C_result);
		double C_g;
	    XL_readNumCell(XL_g,C_g,	        " ARM_ERR: g numeric expected",C_result);
		double C_g0;
	    XL_readNumCell(XL_g0,C_g0,	        " ARM_ERR: g0 numeric expected",C_result);
		double C_AveragingPeriodNb;
	    XL_readNumCell(XL_AveragingPeriodNb,C_AveragingPeriodNb,	        " ARM_ERR: AveragingPeriodNb numeric expected",C_result);
	 	double C_AsianReset;
	    XL_readNumCell(XL_AsianReset,C_AsianReset,	        " ARM_ERR: AsianReset numeric expected",C_result);
		double C_Cashspread;
	    XL_readNumCell(XL_Cashspread,C_Cashspread,	        " ARM_ERR: Cashspread numeric expected",C_result);
	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);
		double C_n;
	    XL_readNumCell(XL_n,C_n,	        " ARM_ERR: n numeric expected",C_result);
		double C_nz;
	    XL_readNumCell(XL_nz,C_nz,	        " ARM_ERR: nz numeric expected",C_result);
		double C_nh;
	    XL_readNumCell(XL_nh,C_nh,	        " ARM_ERR: nz numeric expected",C_result);



	    long retCode = ARMLOCAL_Mepi_VanillaOption_STOBS(
			C_date,
			C_P0,
		    C_K,
		    C_T,
		    C_Zero_Curve_id,
			C_b,
			C_YearlyFees,
			C_Cashspread,
			C_SABR_Model_id,
			C_minExp,
			C_maxExp,
			C_riskFac,
			C_g0,
			C_g,
			C_AveragingPeriodNb,
			C_AsianReset,
		    C_callOrPut,
			C_n,
			C_nz,
			C_nh,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Mepi_VanillaOption_STOBS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_STOBS_delta(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nz,
	LPXLOPER XL_nh)
{
	ADD_LOG("Local_Mepi_VanillaOption_STOBS_delta");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs	
		double C_date;
	    XL_readNumCell(XL_date,C_date,		" ARM_ERR: date numeric expected",C_result);	
	    double C_P0;
	    XL_readNumCell(XL_P0,C_P0,		" ARM_ERR: P0 numeric expected",C_result);	
	    double C_K;
	    XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	    double C_T;
	    XL_readNumCell(XL_T,C_T,	            " ARM_ERR: T numeric expected",C_result);	
		

		long C_Zero_Curve_id;		
		XL_GETOBJID( XL_ZeroCurve, C_Zero_Curve_id,	" ARM_ERR: Zero Curve: Object expected",		C_result);
	
		double C_b;
		XL_readNumCell(XL_b,C_b,	        " ARM_ERR: b numeric expected",C_result);
		double C_YearlyFees;
		XL_readNumCell(XL_YearlyFees,C_YearlyFees,	        " ARM_ERR: YearlyFees numeric expected",C_result);
		
		
		long C_SABR_Model_id;		
		XL_GETOBJID( XL_SABR, C_SABR_Model_id,	" ARM_ERR: SABR Model: Object expected",		C_result);
		
	
		double C_minExp;
		XL_readNumCell(XL_minExp,C_minExp,	        " ARM_ERR: minExp numeric expected",C_result);
		double C_maxExp;
		XL_readNumCell(XL_maxExp,C_maxExp,	        " ARM_ERR: maxExp numeric expected",C_result);
		double C_riskFac;
		XL_readNumCell(XL_riskFac,C_riskFac,	        " ARM_ERR: riskFac numeric expected",C_result);
		double C_g;
	    XL_readNumCell(XL_g,C_g,	        " ARM_ERR: g numeric expected",C_result);
		double C_g0;
	    XL_readNumCell(XL_g0,C_g0,	        " ARM_ERR: g0 numeric expected",C_result);
		double C_AveragingPeriodNb;
	    XL_readNumCell(XL_AveragingPeriodNb,C_AveragingPeriodNb,	        " ARM_ERR: AveragingPeriodNb numeric expected",C_result);
	 	double C_AsianReset;
	    XL_readNumCell(XL_AsianReset,C_AsianReset,	        " ARM_ERR: AsianReset numeric expected",C_result);
		double C_Cashspread;
	    XL_readNumCell(XL_Cashspread,C_Cashspread,	        " ARM_ERR: Cashspread numeric expected",C_result);
	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);
		double C_n;
	    XL_readNumCell(XL_n,C_n,	        " ARM_ERR: n numeric expected",C_result);
		double C_nz;
	    XL_readNumCell(XL_nz,C_nz,	        " ARM_ERR: nz numeric expected",C_result);
		double C_nh;
	    XL_readNumCell(XL_nh,C_nh,	        " ARM_ERR: nz numeric expected",C_result);



	    long retCode = ARMLOCAL_Mepi_VanillaOption_STOBS_delta(
			C_date,
			C_P0,
		    C_K,
		    C_T,
		    C_Zero_Curve_id,
			C_b,
			C_YearlyFees,
			C_Cashspread,
			C_SABR_Model_id,
			C_minExp,
			C_maxExp,
			C_riskFac,
			C_g0,
			C_g,
			C_AveragingPeriodNb,
			C_AsianReset,
		    C_callOrPut,
			C_n,
			C_nz,
			C_nh,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Mepi_VanillaOption_STOBS_delta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_SABR(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_SpotName,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nbMC)
{
	ADD_LOG("Local_Mepi_VanillaOption_SABR");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs	
		double C_date;
	    XL_readNumCell(XL_date,C_date,		" ARM_ERR: date numeric expected",C_result);	
	    double C_P0;
	    XL_readNumCell(XL_P0,C_P0,		" ARM_ERR: P0 numeric expected",C_result);	
	    double C_K;
	    XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	    double C_T;
	    XL_readNumCell(XL_T,C_T,	            " ARM_ERR: T numeric expected",C_result);	
		

		long C_Zero_Curve_id;		
		XL_GETOBJID( XL_ZeroCurve, C_Zero_Curve_id,	" ARM_ERR: Zero Curve: Object expected",		C_result);
		double C_b;
		XL_readNumCell(XL_b,C_b,	        " ARM_ERR: b numeric expected",C_result);
		double C_YearlyFees;
		XL_readNumCell(XL_YearlyFees,C_YearlyFees,	        " ARM_ERR: YearlyFees numeric expected",C_result);
		
		
		long C_SABR_Model_id;		
		XL_GETOBJID( XL_SABR, C_SABR_Model_id,	" ARM_ERR: SABR Model: Object expected",		C_result);
		CCString C_SpotName;
		XL_readStrCell( XL_SpotName, C_SpotName, " ARM_ERR: SpotName: Object expected",C_result);
	
		double C_minExp;
		XL_readNumCell(XL_minExp,C_minExp,	        " ARM_ERR: minExp numeric expected",C_result);
		double C_maxExp;
		XL_readNumCell(XL_maxExp,C_maxExp,	        " ARM_ERR: maxExp numeric expected",C_result);
		double C_riskFac;
		XL_readNumCell(XL_riskFac,C_riskFac,	        " ARM_ERR: riskFac numeric expected",C_result);
		double C_g;
	    XL_readNumCell(XL_g,C_g,	        " ARM_ERR: g numeric expected",C_result);
		double C_g0;
	    XL_readNumCell(XL_g0,C_g0,	        " ARM_ERR: g0 numeric expected",C_result);
		double C_AveragingPeriodNb;
	    XL_readNumCell(XL_AveragingPeriodNb,C_AveragingPeriodNb,	        " ARM_ERR: AveragingPeriodNb numeric expected",C_result);
	 	double C_AsianReset;
	    XL_readNumCell(XL_AsianReset,C_AsianReset,	        " ARM_ERR: AsianReset numeric expected",C_result);
		double C_Cashspread;
	    XL_readNumCell(XL_Cashspread,C_Cashspread,	        " ARM_ERR: Cashspread numeric expected",C_result);
	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);
		double C_n;
	    XL_readNumCell(XL_n,C_n,	        " ARM_ERR: n numeric expected",C_result);
		double C_nbMC;
	    XL_readNumCell(XL_nbMC,C_nbMC,	        " ARM_ERR: nbMC numeric expected",C_result);

		double C_nz=0.0;
	    long retCode = ARMLOCAL_Mepi_VanillaOption_SABR(
			C_date,
			C_P0,
		    C_K,
		    C_T,
		    C_Zero_Curve_id,
			C_b,
			C_YearlyFees,
			C_Cashspread,
			C_SABR_Model_id,
			C_SpotName,
			C_minExp,
			C_maxExp,
			C_riskFac,
			C_g0,
			C_g,
			C_AveragingPeriodNb,
			C_AsianReset,
		    C_callOrPut,
			C_n,
			C_nz,
			C_nbMC,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Mepi_VanillaOption_SABR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_SABR_delta(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_SpotName,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nbMC)
{
	ADD_LOG("Local_Mepi_VanillaOption_SABR_delta");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
	    static char* reason = "";	

	    /// gets the inputs	
		double C_date;
	    XL_readNumCell(XL_date,C_date,		" ARM_ERR: date numeric expected",C_result);	
	    double C_P0;
	    XL_readNumCell(XL_P0,C_P0,		" ARM_ERR: P0 numeric expected",C_result);	
	    double C_K;
	    XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	    double C_T;
	    XL_readNumCell(XL_T,C_T,	            " ARM_ERR: T numeric expected",C_result);	
		

		long C_Zero_Curve_id;		
		XL_GETOBJID( XL_ZeroCurve, C_Zero_Curve_id,	" ARM_ERR: Zero Curve: Object expected",		C_result);
		double C_b;
		XL_readNumCell(XL_b,C_b,	        " ARM_ERR: b numeric expected",C_result);
		double C_YearlyFees;
		XL_readNumCell(XL_YearlyFees,C_YearlyFees,	        " ARM_ERR: YearlyFees numeric expected",C_result);
		
		
		long C_SABR_Model_id;		
		XL_GETOBJID( XL_SABR, C_SABR_Model_id,	" ARM_ERR: SABR Model: Object expected",		C_result);
		CCString C_SpotName;
		XL_readStrCell( XL_SpotName, C_SpotName, " ARM_ERR: SpotName: Object expected",C_result);
	
	
		double C_minExp;
		XL_readNumCell(XL_minExp,C_minExp,	        " ARM_ERR: minExp numeric expected",C_result);
		double C_maxExp;
		XL_readNumCell(XL_maxExp,C_maxExp,	        " ARM_ERR: maxExp numeric expected",C_result);
		double C_riskFac;
		XL_readNumCell(XL_riskFac,C_riskFac,	        " ARM_ERR: riskFac numeric expected",C_result);
		double C_g;
	    XL_readNumCell(XL_g,C_g,	        " ARM_ERR: g numeric expected",C_result);
		double C_g0;
	    XL_readNumCell(XL_g0,C_g0,	        " ARM_ERR: g0 numeric expected",C_result);
		double C_AveragingPeriodNb;
	    XL_readNumCell(XL_AveragingPeriodNb,C_AveragingPeriodNb,	        " ARM_ERR: AveragingPeriodNb numeric expected",C_result);
	 	double C_AsianReset;
	    XL_readNumCell(XL_AsianReset,C_AsianReset,	        " ARM_ERR: AsianReset numeric expected",C_result);
		double C_Cashspread;
	    XL_readNumCell(XL_Cashspread,C_Cashspread,	        " ARM_ERR: Cashspread numeric expected",C_result);	
	    double C_callOrPut;
        XL_GETCONVCALLORPUT(XL_CallOrPut,C_callOrPut,   " ARM_ERR: call or put string expected",C_result);
		double C_n;
	    XL_readNumCell(XL_n,C_n,	        " ARM_ERR: n numeric expected",C_result);
			double C_nbMC;
	    XL_readNumCell(XL_nbMC,C_nbMC,	        " ARM_ERR: nbMC numeric expected",C_result);
		double C_nz=0.0;
	    long retCode = ARMLOCAL_Mepi_VanillaOption_SABR_delta(
			C_date,
			C_P0,
		    C_K,
		    C_T,
		    C_Zero_Curve_id,
			C_b,
			C_YearlyFees,
			C_Cashspread,
			C_SABR_Model_id,
			C_SpotName,
			C_minExp,
			C_maxExp,
			C_riskFac,
			C_g0,
			C_g,
			C_AveragingPeriodNb,
			C_AsianReset,
		    C_callOrPut,
			C_n,
			C_nz,
			C_nbMC,
		    C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Mepi_VanillaOption_SABR_delta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Fund VanillaOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////




__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_averaging,	
	LPXLOPER XL_reset,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_StochasticVol_LN_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_averaging;
	XL_readNumCell(XL_averaging,C_averaging,	" ARM_ERR: averaging numeric expected",C_result);
	double C_reset;
	XL_readNumCell(XL_reset,C_reset,	" ARM_ERR: reset numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_LN_VanillaOption(
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_averaging,
		C_reset,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_averaging,	
	LPXLOPER XL_reset,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_StochasticVol_LN_VanillaOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_averaging;
	XL_readNumCell(XL_averaging,C_averaging,	" ARM_ERR: averaging numeric expected",C_result);
	double C_reset;
	XL_readNumCell(XL_reset,C_reset,	" ARM_ERR: reset numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_LN_VanillaOption_Der(
		C_i,
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_averaging,
		C_reset,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_averaging,	
	LPXLOPER XL_reset,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_StochasticVol_LN_VanillaOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_averaging;
	XL_readNumCell(XL_averaging,C_averaging,	" ARM_ERR: averaging numeric expected",C_result);
	double C_reset;
	XL_readNumCell(XL_reset,C_reset,	" ARM_ERR: reset numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_LN_VanillaOption_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_averaging,
		C_reset,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_Ari_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_averaging,	
	LPXLOPER XL_reset,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_StochasticVol_LN_Ari_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_averaging;
	XL_readNumCell(XL_averaging,C_averaging,	" ARM_ERR: averaging numeric expected",C_result);
	double C_reset;
	XL_readNumCell(XL_reset,C_reset,	" ARM_ERR: reset numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_LN_Ari_VanillaOption(
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_averaging,
		C_reset,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_Ari_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_averaging,	
	LPXLOPER XL_reset,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_StochasticVol_LN_Ari_VanillaOption_Der");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_averaging;
	XL_readNumCell(XL_averaging,C_averaging,	" ARM_ERR: averaging numeric expected",C_result);
	double C_reset;
	XL_readNumCell(XL_reset,C_reset,	" ARM_ERR: reset numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_LN_Ari_VanillaOption_Der(
		C_i,
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_averaging,
		C_reset,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_Ari_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_averaging,	
	LPXLOPER XL_reset,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_StochasticVol_LN_Ari_VanillaOption_Der2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_i;
	XL_readNumCell(XL_i,C_i,		" ARM_ERR: i numeric expected",C_result);
	double C_j;
	XL_readNumCell(XL_j,C_j,		" ARM_ERR: j numeric expected",C_result);
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: sig numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_averaging;
	XL_readNumCell(XL_averaging,C_averaging,	" ARM_ERR: averaging numeric expected",C_result);
	double C_reset;
	XL_readNumCell(XL_reset,C_reset,	" ARM_ERR: reset numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_LN_Ari_VanillaOption_Der2(
		C_i,
		C_j,
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_averaging,
		C_reset,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///			OptimizeSkewVector
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_OptimizeSkewVectorCommon(
	LPXLOPER XL_C_tagetVect,
	LPXLOPER XL_C_weights,
    LPXLOPER XL_C_presicions,
	LPXLOPER XL_C_InitVector,
	LPXLOPER XL_C_LBoundVector,
	LPXLOPER XL_C_UBoundVector,
	LPXLOPER XL_C_algo,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		VECTOR<double> C_tagetVect;
		XL_readNumVector(XL_C_tagetVect,C_tagetVect," ARM_ERR: target Vector: array of numeric expected",C_result);

		VECTOR<double> C_weights;
		XL_readNumVector(XL_C_weights,C_weights," ARM_ERR: target Vector: array of numeric expected",C_result);

        VECTOR<double> C_presicions;
		XL_readNumVector(XL_C_presicions,C_presicions," ARM_ERR: precisions Vector: array of numeric expected",C_result);

		VECTOR<double> C_InitVector;
		XL_readNumVector(XL_C_InitVector,C_InitVector," ARM_ERR: target Vector: array of numeric expected",C_result);

		VECTOR<double> C_LBoundVector;
		XL_readNumVector(XL_C_LBoundVector,C_LBoundVector," ARM_ERR: target Vector: array of numeric expected",C_result);

		VECTOR<double> C_UBoundVector;
		XL_readNumVector(XL_C_UBoundVector,C_UBoundVector," ARM_ERR: target Vector: array of numeric expected",C_result);
		
		/// Vector to store the result
		vector<double> C_DataResult;

		double C_algoDble;
		XL_readNumCell( XL_C_algo, C_algoDble, " ARM_ERR: new abscisse: numeric expected",	C_result);
		long C_algo = (long) C_algoDble;
		/// a function with a context

		ARMLOCAL_OptimizeSkewVector(C_tagetVect,C_weights,C_presicions,C_InitVector,C_LBoundVector,C_UBoundVector,C_DataResult,C_algo,C_result);
		/// add these additional lines 
		/// to display blank lines
		const int additionalLinesNb = 100;
		bool fillWithBlank = true;
		FreeCurCellContent ();
		XL_writeNumVectorWithOptions( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PricingModel_GetModelMapCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_OptimizeSkewVector(
		LPXLOPER XL_C_tagetVect,
		LPXLOPER XL_C_weights,
        LPXLOPER XL_C_presicions,
		LPXLOPER XL_C_InitVector,
		LPXLOPER XL_C_LBoundVector,
		LPXLOPER XL_C_UBoundVector,
		LPXLOPER XL_C_algo)
{
	bool PersistentInXL = true;
	return Local_OptimizeSkewVectorCommon(
		XL_C_tagetVect,
		XL_C_weights,
        XL_C_presicions,
		XL_C_InitVector,
		XL_C_LBoundVector,
		XL_C_UBoundVector,
		XL_C_algo,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_OptimizeSkewVector(
		LPXLOPER XL_C_tagetVect,
		LPXLOPER XL_C_weights,
        LPXLOPER XL_C_presicions,
		LPXLOPER XL_C_InitVector,
		LPXLOPER XL_C_LBoundVector,
		LPXLOPER XL_C_UBoundVector,
		LPXLOPER XL_C_algo)
{
	ADD_LOG("Local_PXL_OptimizeSkewVector");
	bool PersistentInXL = false;
	return Local_OptimizeSkewVectorCommon(
		XL_C_tagetVect,
		XL_C_weights,
        XL_C_presicions,
		XL_C_InitVector,
		XL_C_LBoundVector,
		XL_C_UBoundVector,
		XL_C_algo,
		PersistentInXL );
}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Stochastic Normal volatility Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_N_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	)



{
	ADD_LOG("Local_StochasticVol_N_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs	
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);

	double C_r;
	XL_readNumCell(XL_r,C_r, " ARM_ERR: r numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,	" ARM_ERR: stddev numeric expected",C_result);	
	double C_VolDrift;
	XL_readNumCell(XL_VolDrift,C_VolDrift,	" ARM_ERR: VolDrift numeric expected",C_result);
	double C_VolVol;
	XL_readNumCell(XL_VolVol,C_VolVol,	" ARM_ERR: VolVol numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);	
	double C_nb;
	double n_default;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_StochasticVol_N_VanillaOption(
		C_f,
		C_K,
		C_T,
		C_r,
		C_sig,
		C_VolDrift,
		C_VolVol,
		C_callput,
		C_nb,
		C_result
	);
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       GLambda Distribution
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_From_SABR_Calibrate(
	
	LPXLOPER XL_forward,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_T,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_Scope,
	LPXLOPER XL_IniL1,
	LPXLOPER XL_IniL2,
	LPXLOPER XL_IniL3,
	LPXLOPER XL_IniL4,
	LPXLOPER XL_IniL5,
	LPXLOPER XL_IniL6,
	LPXLOPER XL_Algo
	)
{
	ADD_LOG("Local_GLambda_From_SABR_Calibrate");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		double C_forward;
		XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward Value numeric expected",C_result);
		double C_alpha;
		XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha Value numeric expected",C_result);
		double C_beta;
		XL_readNumCell(XL_beta,C_beta,		" ARM_ERR: beta Value numeric expected",C_result);
		double C_rho;
		XL_readNumCell(XL_rho,C_rho,		" ARM_ERR: rho Value numeric expected",C_result);
		double C_nu;
		XL_readNumCell(XL_nu,C_nu,		" ARM_ERR: nu Value numeric expected",C_result);
		double C_T;
		XL_readNumCell(XL_T,C_T,		" ARM_ERR: T Value numeric expected",C_result);
		double C_SabrType;
		XL_GETCONVSABRFLAG(XL_SabrType,C_SabrType," ARM_ERR: SabrType flag string expected",C_result);
		double C_Scope;
		XL_readNumCell(XL_Scope,C_Scope,		" ARM_ERR: Scope Value numeric expected",C_result);
		double C_IniL1;
		XL_readNumCell(XL_IniL1,C_IniL1,		" ARM_ERR: IniL1 Value numeric expected",C_result);
		double C_IniL2;
		XL_readNumCell(XL_IniL2,C_IniL2,		" ARM_ERR: IniL2 Value numeric expected",C_result);
		double C_IniL3;
		XL_readNumCell(XL_IniL3,C_IniL3,		" ARM_ERR: IniL3 Value numeric expected",C_result);
		double C_IniL4;
		XL_readNumCell(XL_IniL4,C_IniL4,		" ARM_ERR: IniL4 Value numeric expected",C_result);
		double C_IniL5;
		XL_readNumCell(XL_IniL5,C_IniL5,		" ARM_ERR: IniL5 Value numeric expected",C_result);
		double C_IniL6;
		XL_readNumCell(XL_IniL6,C_IniL6,		" ARM_ERR: IniL6 Value numeric expected",C_result);
		double C_Algo;
		XL_GETCONVOPTIMIZATIONFLAG(XL_Algo,C_Algo," ARM_ERR: Algo flag string expected",C_result);
	
		/// call the function
		long retCode=ARMLOCAL_GLambda_From_SABR_Calibrate(
			C_forward,
			C_alpha,
			C_beta,
			C_rho,
			C_nu,
			C_T,
			C_SabrType,
			C_Scope,
			C_IniL1,
			C_IniL2,
			C_IniL3,
			C_IniL4,
			C_IniL5,
			C_IniL6,
			C_Algo,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 7;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GLambda_From_SABR_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

	

/////////////////////////////////////////////////////////////////////////////////////
///
///		GLambda + Student law
///
/////////////////////////////////////////////////////////////////////////////////////




__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Student_GLambda_Power_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	

	double C_l1b;
	XL_readNumCell(XL_l1b,C_l1b,		" ARM_ERR: l1b numeric expected",C_result);	
	double C_l2b;
	XL_readNumCell(XL_l2b,C_l2b,		" ARM_ERR: l2b numeric expected",C_result);	
	double C_l3b;
	XL_readNumCell(XL_l3b,C_l3b,		" ARM_ERR: l3b numeric expected",C_result);	
	double C_l4b;
	XL_readNumCell(XL_l4b,C_l4b,		" ARM_ERR: l4b numeric expected",C_result);	
	double C_l5b;
	XL_readNumCell(XL_l5b,C_l5b,		" ARM_ERR: l5b numeric expected",C_result);	
	double C_l6b;
	XL_readNumCell(XL_l6b,C_l6b,		" ARM_ERR: l6b numeric expected",C_result);	

	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_copula_degre;
	XL_readNumCell(XL_copula_degre,C_copula_degre,		" ARM_ERR: copula_degre numeric expected",C_result);	

	double C_a10;
	C_a10=1.0;  /// due to limitation of the excel interface, this coeficient is put to 1

	double C_b10;
	XL_readNumCell(XL_b10,C_b10,		" ARM_ERR: b10 numeric expected",C_result);
	double C_k10;
	XL_readNumCell(XL_k10,C_k10,		" ARM_ERR: k10 numeric expected",C_result);
	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default=40;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Student_GLambda_Power_SpreadOption(
		C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_l1b,
		C_l2b,
		C_l3b,
		C_l4b,
		C_l5b,
		C_l6b,
		C_copula_corr,
		C_copula_degre,
		C_a10,
		C_b10,
		C_k10,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_GLambda_Power_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the gamma function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Lambert_Function(
	LPXLOPER XL_x
	)
	{
		ADD_LOG("Local_Lambert_Function");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Lambert_Function(
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Lambert_Function" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Asymptotic Time Value Formula
///
////////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholesTimeValue(
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	)
	{
		ADD_LOG("Local_CF_BlackSholesTimeValue");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_totalvolatility;
	XL_readNumCell(XL_totalvolatility,C_totalvolatility,		" ARM_ERR: totalvolatility numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholesTimeValue(
		C_forward,
		C_totalvolatility,
		C_bondprice,
		C_strike,
		C_callput,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholesTimeValue_ImplicitVol(
	LPXLOPER XL_forward,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_optprice
	)
	{
		ADD_LOG("Local_CF_BlackSholesTimeValue_ImplicitVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,		" ARM_ERR: forward numeric expected",C_result);	
	double C_bondprice;
	XL_readNumCell(XL_bondprice,C_bondprice,	" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,	" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_optprice;
	XL_readNumCell(XL_optprice,C_optprice,		" ARM_ERR: forward numeric expected",C_result);	

	
	/// call the function
	long retCode=ARMLOCAL_CF_BlackSholesTimeValue_ImplicitVol(
		C_forward,
		C_bondprice,
		C_strike,
		C_callput,
		C_optprice,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CF_BlackSholes_ImplicitVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedLogNormal_Quantile(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_alpha
	)
	{
		ADD_LOG("Local_ShiftedLogNormal_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,	" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: T numeric expected",C_result);	

	double C_sigma;
	XL_readNumCell(XL_sigma,C_sigma,		" ARM_ERR: sigma numeric expected",C_result);		
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	

	
	/// call the function
	long retCode=ARMLOCAL_CF_ShiftedLogNormal_Quantile(
		C_f,
		C_k,
		C_t,
		C_sigma,
		C_alpha,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ShiftedLogNormal_Quantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedLogNormal_Distribution(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_alpha
	)
	{
		ADD_LOG("Local_ShiftedLogNormal_Distribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,	" ARM_ERR: k numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,	" ARM_ERR: T numeric expected",C_result);	

	double C_sigma;
	XL_readNumCell(XL_sigma,C_sigma,		" ARM_ERR: sigma numeric expected",C_result);		
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	

	
	/// call the function
	long retCode=ARMLOCAL_CF_ShiftedLogNormal_Distribution(
		C_f,
		C_k,
		C_t,
		C_sigma,
		C_alpha,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ShiftedLogNormal_Distribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Quantile(
	LPXLOPER XL_f,
	LPXLOPER XL_strike,
	LPXLOPER XL_T,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
		double C_strike;
	XL_readNumCell(XL_strike,C_strike,		" ARM_ERR: strike numeric expected",C_result);

		double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	

	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,	" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	

	double C_rho;
	XL_readNumCell(XL_rho,C_rho,		" ARM_ERR: rho numeric expected",C_result);		
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,		" ARM_ERR: nu numeric expected",C_result);	
	double C_SabrType;
	XL_GETCONVSABRFLAG(XL_SabrType,C_SabrType," ARM_ERR: SabrType flag string expected",C_result);
	
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	
	/// call the function
	long retCode=ARMLOCAL_CF_SABR_Quantile(
		C_f,
		C_strike,
		C_T,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_SabrType,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARMLOCAL_CF_SABR_Quantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Distribution(
	LPXLOPER XL_f,
	LPXLOPER XL_strike,
	LPXLOPER XL_T,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	)
	{
		ADD_LOG("Local_SABR_Distribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f;
	XL_readNumCell(XL_f,C_f,		" ARM_ERR: f numeric expected",C_result);
		double C_strike;
	XL_readNumCell(XL_strike,C_strike,		" ARM_ERR: strike numeric expected",C_result);

		double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	

	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,	" ARM_ERR: alpha numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	

	double C_rho;
	XL_readNumCell(XL_rho,C_rho,		" ARM_ERR: rho numeric expected",C_result);		
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,		" ARM_ERR: nu numeric expected",C_result);	
	double C_SabrType;
	XL_GETCONVSABRFLAG(XL_SabrType,C_SabrType," ARM_ERR: SabrType flag string expected",C_result);
	
	double C_nbsteps;
	double n_default=120;
    XL_readNumCellWD(XL_nbsteps, C_nbsteps, n_default, " ARM_ERR: nbsteps: numeric expected",C_result);
	double WD_Alpha_Exp=0.01;
	double C_Alpha_Exp;
    XL_readNumCellWD(XL_Alpha_Exp, C_Alpha_Exp, WD_Alpha_Exp, " ARM_ERR: WD_Alpha_Exp: numeric expected",C_result);
	double WD_Alpha_Tanh=0.03;
	double C_Alpha_Tanh;
    XL_readNumCellWD(XL_Alpha_Tanh, C_Alpha_Tanh, WD_Alpha_Tanh, " ARM_ERR: WD_Alpha_Tanh: numeric expected",C_result);
	double WD_Kb_Tanh=1.5;
	double C_Kb_Tanh;
    XL_readNumCellWD(XL_Kb_Tanh, C_Kb_Tanh, WD_Kb_Tanh, " ARM_ERR: Kb_Tanh: numeric expected",C_result);

	
	/// call the function
	long retCode=ARMLOCAL_CF_SABR_Distribution(
		C_f,
		C_strike,
		C_T,
		C_alpha,
		C_beta,
		C_rho,
		C_nu,
		C_SabrType,
		C_nbsteps,
		C_Alpha_Exp,
		C_Alpha_Tanh,
		C_Kb_Tanh,
		C_result
		);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARMLOCAL_CF_SABR_Distribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/////////////////////////////////////////////////////////////////////////////////////
///
/// Digital Spread Option Student on Glambda underlyings
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_Digital_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Student_GLambda_Power_Digital_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	

	double C_l1b;
	XL_readNumCell(XL_l1b,C_l1b,		" ARM_ERR: l1b numeric expected",C_result);	
	double C_l2b;
	XL_readNumCell(XL_l2b,C_l2b,		" ARM_ERR: l2b numeric expected",C_result);	
	double C_l3b;
	XL_readNumCell(XL_l3b,C_l3b,		" ARM_ERR: l3b numeric expected",C_result);	
	double C_l4b;
	XL_readNumCell(XL_l4b,C_l4b,		" ARM_ERR: l4b numeric expected",C_result);	
	double C_l5b;
	XL_readNumCell(XL_l5b,C_l5b,		" ARM_ERR: l5b numeric expected",C_result);	
	double C_l6b;
	XL_readNumCell(XL_l6b,C_l6b,		" ARM_ERR: l6b numeric expected",C_result);	

	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_copula_degre;
	XL_readNumCell(XL_copula_degre,C_copula_degre,		" ARM_ERR: copula_degre numeric expected",C_result);	

	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Student_GLambda_Power_Digital_SpreadOption(
			C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_l1b,
		C_l2b,
		C_l3b,
		C_l4b,
		C_l5b,
		C_l6b,
		C_copula_corr,
		C_copula_degre,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_GLambda_Power_Digital_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_Index1Digital_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Student_GLambda_Power_Index1Digital_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	

	double C_l1b;
	XL_readNumCell(XL_l1b,C_l1b,		" ARM_ERR: l1b numeric expected",C_result);	
	double C_l2b;
	XL_readNumCell(XL_l2b,C_l2b,		" ARM_ERR: l2b numeric expected",C_result);	
	double C_l3b;
	XL_readNumCell(XL_l3b,C_l3b,		" ARM_ERR: l3b numeric expected",C_result);	
	double C_l4b;
	XL_readNumCell(XL_l4b,C_l4b,		" ARM_ERR: l4b numeric expected",C_result);	
	double C_l5b;
	XL_readNumCell(XL_l5b,C_l5b,		" ARM_ERR: l5b numeric expected",C_result);	
	double C_l6b;
	XL_readNumCell(XL_l6b,C_l6b,		" ARM_ERR: l6b numeric expected",C_result);	

	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_copula_degre;
	XL_readNumCell(XL_copula_degre,C_copula_degre,		" ARM_ERR: copula_degre numeric expected",C_result);	

	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Student_GLambda_Power_Index1Digital_SpreadOption(
			C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_l1b,
		C_l2b,
		C_l3b,
		C_l4b,
		C_l5b,
		C_l6b,
		C_copula_corr,
		C_copula_degre,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_GLambda_Power_Index1Digital_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_Index2Digital_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Student_GLambda_Power_Index2Digital_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	

	double C_l1b;
	XL_readNumCell(XL_l1b,C_l1b,		" ARM_ERR: l1b numeric expected",C_result);	
	double C_l2b;
	XL_readNumCell(XL_l2b,C_l2b,		" ARM_ERR: l2b numeric expected",C_result);	
	double C_l3b;
	XL_readNumCell(XL_l3b,C_l3b,		" ARM_ERR: l3b numeric expected",C_result);	
	double C_l4b;
	XL_readNumCell(XL_l4b,C_l4b,		" ARM_ERR: l4b numeric expected",C_result);	
	double C_l5b;
	XL_readNumCell(XL_l5b,C_l5b,		" ARM_ERR: l5b numeric expected",C_result);	
	double C_l6b;
	XL_readNumCell(XL_l6b,C_l6b,		" ARM_ERR: l6b numeric expected",C_result);	

	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_copula_degre;
	XL_readNumCell(XL_copula_degre,C_copula_degre,		" ARM_ERR: copula_degre numeric expected",C_result);	

	double C_a20;
	XL_readNumCell(XL_a20,C_a20,		" ARM_ERR: a20 numeric expected",C_result);
	double C_b20;
	XL_readNumCell(XL_b20,C_b20,		" ARM_ERR: b20 numeric expected",C_result);
	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Student_GLambda_Power_Index2Digital_SpreadOption(
			C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_l1b,
		C_l2b,
		C_l3b,
		C_l4b,
		C_l5b,
		C_l6b,
		C_copula_corr,
		C_copula_degre,
		C_a20,
		C_b20,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_GLambda_Power_Index2Digital_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_Student_GLambda_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	

	double C_l1b;
	XL_readNumCell(XL_l1b,C_l1b,		" ARM_ERR: l1b numeric expected",C_result);	
	double C_l2b;
	XL_readNumCell(XL_l2b,C_l2b,		" ARM_ERR: l2b numeric expected",C_result);	
	double C_l3b;
	XL_readNumCell(XL_l3b,C_l3b,		" ARM_ERR: l3b numeric expected",C_result);	
	double C_l4b;
	XL_readNumCell(XL_l4b,C_l4b,		" ARM_ERR: l4b numeric expected",C_result);	
	double C_l5b;
	XL_readNumCell(XL_l5b,C_l5b,		" ARM_ERR: l5b numeric expected",C_result);	
	double C_l6b;
	XL_readNumCell(XL_l6b,C_l6b,		" ARM_ERR: l6b numeric expected",C_result);	

	double C_copula_corr;
	XL_readNumCell(XL_copula_corr,C_copula_corr,		" ARM_ERR: copula_corr numeric expected",C_result);	
	double C_copula_degre;
	XL_readNumCell(XL_copula_degre,C_copula_degre,		" ARM_ERR: copula_degre numeric expected",C_result);	

	double C_k20;
	XL_readNumCell(XL_k20,C_k20,		" ARM_ERR: k20 numeric expected",C_result);
	double C_n;
	double n_default;
    XL_readNumCellWD(XL_n, C_n, n_default, " ARM_ERR: n: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Student_GLambda_SpreadOption(
			C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_l1b,
		C_l2b,
		C_l3b,
		C_l4b,
		C_l5b,
		C_l6b,
		C_copula_corr,
		C_copula_degre,
		C_k20,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_GLambda_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_Distribution(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_GLambda_Distribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	


	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_GLambda_Distribution(
			C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GLambda_Distribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_Quantile(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_GLambda_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_l1a;
	XL_readNumCell(XL_l1a,C_l1a,		" ARM_ERR: l1a numeric expected",C_result);	
	double C_l2a;
	XL_readNumCell(XL_l2a,C_l2a,		" ARM_ERR: l2a numeric expected",C_result);	
	double C_l3a;
	XL_readNumCell(XL_l3a,C_l3a,		" ARM_ERR: l3a numeric expected",C_result);	
	double C_l4a;
	XL_readNumCell(XL_l4a,C_l4a,		" ARM_ERR: l4a numeric expected",C_result);	
	double C_l5a;
	XL_readNumCell(XL_l5a,C_l5a,		" ARM_ERR: l5a numeric expected",C_result);	
	double C_l6a;
	XL_readNumCell(XL_l6a,C_l6a,		" ARM_ERR: l6a numeric expected",C_result);	


	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_GLambda_Quantile(
			C_l1a,
		C_l2a,
		C_l3a,
		C_l4a,
		C_l5a,
		C_l6a,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GLambda_Quantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Student_Quantile(
	LPXLOPER XL_degre,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Student_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_degre;
	XL_readNumCell(XL_degre,C_degre,		" ARM_ERR: degre numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_Student_Quantile(
		C_degre,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_Quantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Student_Distribution(
	LPXLOPER XL_degre,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Student_Distribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_degre;
	XL_readNumCell(XL_degre,C_degre,		" ARM_ERR: degre numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_Student_Distribution(
		C_degre,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Student_Distribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ImcompleteBeta_Inverse(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_ImcompleteBeta_Inverse");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_ImcompleteBeta_Inverse(
		C_a,
		C_b,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ImcompleteBeta_Inverse" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Student_QIntegral(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Student_QIntegral");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_Student_QIntegral(
		C_a,
		C_b,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ImcompleteBeta_Inverse" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Normal_ImpliedVol(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x,
	LPXLOPER XL_CallPut
	)
{
	ADD_LOG("Local_Normal_ImpliedVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_Normal_ImpliedVol(
		C_a,
		C_b,
		C_x,
		C_callput,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_ImpliedVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_Digital_ImpliedVol(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x,
	LPXLOPER XL_CallPut
	)
{
	ADD_LOG("Local_Normal_Digital_ImpliedVol");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);
	double C_callput;
	XL_GETCONVCALLORPUT(XL_CallPut,C_callput," ARM_ERR: call or put string expected",C_result);
	/// call the function
	long retCode=ARMLOCAL_Normal_Digital_ImpliedVol(
		C_a,
		C_b,
		C_x,
		C_callput,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_ImpliedVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric_Whittaker_M(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Hypergeometric_Whittaker_M");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Hypergeometric_Whittaker_M(
		C_a,
		C_b,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Hypergeometric_Whittaker_M" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric_Whittaker_W(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Hypergeometric_Whittaker_W");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_b;
	XL_readNumCell(XL_b,C_b,		" ARM_ERR: b numeric expected",C_result);	

	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Hypergeometric_Whittaker_W(
		C_a,
		C_b,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Hypergeometric_Whittaker_W" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_Y(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Bessel_Y");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Bessel_Y(
		C_a,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bessel_Y" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_I(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Bessel_I");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Bessel_I(
		C_a,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bessel_I" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_J(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Bessel_J");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Bessel_J(
		C_a,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bessel_J" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_K(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	)
{
	ADD_LOG("Local_Bessel_K");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,		" ARM_ERR: a numeric expected",C_result);	
	double C_x;
	XL_readNumCell(XL_x,C_x,		" ARM_ERR: x numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_Bessel_K(
		C_a,
		C_x,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bessel_K" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		Shifted Heston  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_Shifted_Heston_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_kappa,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_Shift,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb1,
	LPXLOPER XL_nb,
	LPXLOPER XL_nbS,
	LPXLOPER XL_nbO,
	LPXLOPER XL_prec
	)
		{
			ADD_LOG("Local_Shifted_Heston_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,		" ARM_ERR: sig numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_kappa;
	XL_readNumCell(XL_kappa,C_kappa,		" ARM_ERR: kappa numeric expected",C_result);	
	double C_theta;
	XL_readNumCell(XL_theta,C_theta,	" ARM_ERR: theta numeric expected",C_result);	
	double C_ksi;
	XL_readNumCell(XL_ksi,C_ksi,	" ARM_ERR: ksi numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_Shift;
	XL_readNumCell(XL_Shift,C_Shift,	" ARM_ERR: Shift numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb1;
	double n1_default;
    XL_readNumCellWD(XL_nb1, C_nb1, n1_default, " ARM_ERR: nb1: numeric expected",C_result);
	double C_nb;
	double n_default=120;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);
	double C_nbS;
	double nS_default=7;
    XL_readNumCellWD(XL_nbS, C_nbS, nS_default, " ARM_ERR: nbS: numeric expected",C_result);
	double C_nbO;
	double nO_default=3;
    XL_readNumCellWD(XL_nbO, C_nbO, nO_default, " ARM_ERR: nbO: numeric expected",C_result);
	double C_prec;
	double np_default=1e-5;
    XL_readNumCellWD(XL_prec, C_prec, np_default, " ARM_ERR: prec: numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_Shifted_Heston_VanillaOption(
		C_F,
		C_K,
		C_sig,
		C_t,
		C_kappa,
		C_theta,
		C_ksi,
		C_rho,
		C_Shift,
		C_callput,
		C_nb1,
		C_nb,
		C_nbS,
		C_nbO,
		C_prec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Shifted_Heston_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR Heston  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Heston_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_kappa,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_beta,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb1,
	LPXLOPER XL_nb,
	LPXLOPER XL_nbS,
	LPXLOPER XL_nbO,
	LPXLOPER XL_prec
	)
		{
			ADD_LOG("Local_SABR_Heston_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F;
	XL_readNumCell(XL_F,C_F,		" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,		" ARM_ERR: sig numeric expected",C_result);	
	double C_t;
	XL_readNumCell(XL_t,C_t,		" ARM_ERR: tex numeric expected",C_result);	
	double C_kappa;
	XL_readNumCell(XL_kappa,C_kappa,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_theta;
	XL_readNumCell(XL_theta,C_theta,	" ARM_ERR: theta numeric expected",C_result);	
	double C_ksi;
	XL_readNumCell(XL_ksi,C_ksi,	" ARM_ERR: ksi numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,	" ARM_ERR: rho numeric expected",C_result);	
	double C_beta;
	XL_readNumCell(XL_beta,C_beta,	" ARM_ERR: beta numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nb1;
	double n1_default;
    XL_readNumCellWD(XL_nb1, C_nb1, n1_default, " ARM_ERR: nb1: numeric expected",C_result);
	double C_nb;
	double n_default=120;
    XL_readNumCellWD(XL_nb, C_nb, n_default, " ARM_ERR: nb: numeric expected",C_result);
	double C_nbS;
	double nS_default=7;
    XL_readNumCellWD(XL_nbS, C_nbS, nS_default, " ARM_ERR: nbS: numeric expected",C_result);
	double C_nbO;
	double nO_default=3;
    XL_readNumCellWD(XL_nbO, C_nbO, nO_default, " ARM_ERR: nbO: numeric expected",C_result);
	double C_prec;
	double np_default=1e-5;
    XL_readNumCellWD(XL_prec, C_prec, np_default, " ARM_ERR: prec: numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_SABR_Heston_VanillaOption(
		C_F,
		C_K,
		C_sig,
		C_t,
		C_kappa,
		C_theta,
		C_ksi,
		C_rho,
		C_beta,
		C_callput,
		C_nb1,
		C_nb,
		C_nbS,
		C_nbO,
		C_prec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Heston_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_CompleteSpreadoption(
	
	LPXLOPER XL_l1a_Vec,
	LPXLOPER XL_l2a_Vec,
	LPXLOPER XL_l3a_Vec,
	LPXLOPER XL_l4a_Vec,
	LPXLOPER XL_l5a_Vec,
	LPXLOPER XL_l6a_Vec,
	LPXLOPER XL_l1b_Vec,
	LPXLOPER XL_l2b_Vec,
	LPXLOPER XL_l3b_Vec,
	LPXLOPER XL_l4b_Vec,
	LPXLOPER XL_l5b_Vec,
	LPXLOPER XL_l6b_Vec,
	LPXLOPER XL_Discount_Vec,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_k,
	LPXLOPER XL_n
	)
{
	ADD_LOG("Local_GLambda_CompleteSpreadoption");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_l1a_Vec;
		XL_readNumVector(XL_l1a_Vec,C_l1a_Vec," ARM_ERR: l1a Vector: array of numeric expected",C_result);
		vector<double> C_l2a_Vec;
		XL_readNumVector(XL_l2a_Vec,C_l2a_Vec," ARM_ERR: l2a Vector: array of numeric expected",C_result);
		vector<double> C_l3a_Vec;
		XL_readNumVector(XL_l3a_Vec,C_l3a_Vec," ARM_ERR: l3a Vector: array of numeric expected",C_result);
		vector<double> C_l4a_Vec;
		XL_readNumVector(XL_l4a_Vec,C_l4a_Vec," ARM_ERR: l4a Vector: array of numeric expected",C_result);
		vector<double> C_l5a_Vec;
		XL_readNumVector(XL_l5a_Vec,C_l5a_Vec," ARM_ERR: l5a Vector: array of numeric expected",C_result);
		vector<double> C_l6a_Vec;
		XL_readNumVector(XL_l6a_Vec,C_l6a_Vec," ARM_ERR: l6a Vector: array of numeric expected",C_result);
		vector<double> C_l1b_Vec;
		XL_readNumVector(XL_l1b_Vec,C_l1b_Vec," ARM_ERR: l1b Vector: array of numeric expected",C_result);
		vector<double> C_l2b_Vec;
		XL_readNumVector(XL_l2b_Vec,C_l2b_Vec," ARM_ERR: l2b Vector: array of numeric expected",C_result);
		vector<double> C_l3b_Vec;
		XL_readNumVector(XL_l3b_Vec,C_l3b_Vec," ARM_ERR: l3b Vector: array of numeric expected",C_result);
		vector<double> C_l4b_Vec;
		XL_readNumVector(XL_l4b_Vec,C_l4b_Vec," ARM_ERR: l4b Vector: array of numeric expected",C_result);
		vector<double> C_l5b_Vec;
		XL_readNumVector(XL_l5b_Vec,C_l5b_Vec," ARM_ERR: l5b Vector: array of numeric expected",C_result);
		vector<double> C_l6b_Vec;
		XL_readNumVector(XL_l6b_Vec,C_l6b_Vec," ARM_ERR: l6b Vector: array of numeric expected",C_result);
		
		vector<double> C_Discount_Vec;
		XL_readNumVector(XL_Discount_Vec,C_Discount_Vec," ARM_ERR: Discount_Vec Vector: array of numeric expected",C_result);
		
		double C_copula_corr;
		XL_readNumCell(XL_copula_corr,C_copula_corr,	" ARM_ERR: copula_corr numeric expected",C_result);	
		
		double C_copula_degre;
		XL_readNumCell(XL_copula_degre,C_copula_degre,	" ARM_ERR: copula_degre numeric expected",C_result);	
		
		double C_k;
		XL_readNumCell(XL_k,C_k,	" ARM_ERR: k numeric expected",C_result);
		
		double np_default=1e-5;
		double C_n;
		XL_readNumCellWD(XL_n, C_n, np_default, " ARM_ERR: n: numeric expected",C_result);

		
		/// call the function
		long retCode=ARMLOCAL_GLambda_CompleteSpreadoption(
			C_l1a_Vec,
			C_l2a_Vec,
			C_l3a_Vec,
			C_l4a_Vec,
			C_l5a_Vec,
			C_l6a_Vec,
			C_l1b_Vec,
			C_l2b_Vec,
			C_l3b_Vec,
			C_l4b_Vec,
			C_l5b_Vec,
			C_l6b_Vec,
			C_Discount_Vec,
			C_copula_corr,
			C_copula_degre,
			C_k,
			C_n,
			C_result);
		
/*		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
			}
			else
			{
			ARM_ERR();
			}
		*/
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GLambda_CompleteSpreadoption" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		JumpDiffusion Mepi  Vanilla option Formula : return price and delta
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_JumpDiffusion_Mepi_Call(
	LPXLOPER XL_P0,
	LPXLOPER XL_f0,
	LPXLOPER XL_T,
	LPXLOPER XL_K,
	LPXLOPER XL_R,
	LPXLOPER XL_Emin,
	LPXLOPER XL_Lmax,
	LPXLOPER XL_gamma0,
	LPXLOPER XL_gamma1,
	LPXLOPER XL_sig,
	LPXLOPER XL_lambda,
	LPXLOPER XL_sigJ,
	LPXLOPER XL_r,
	LPXLOPER XL_s,
	LPXLOPER XL_mu,
	LPXLOPER XL_fees,
	LPXLOPER XL_volDrift,
	LPXLOPER XL_volVol,
	LPXLOPER XL_callput,
	LPXLOPER XL_params
	)
		{
			ADD_LOG("Local_JumpDiffusion_Mepi_Call");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_P0;
	XL_readNumCell(XL_P0,C_P0,				" ARM_ERR: P0 numeric expected",C_result);
	double C_f0;
	XL_readNumCell(XL_f0,C_f0,				" ARM_ERR: f0 numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);	
	double C_R;
	XL_readNumCell(XL_R,C_R,				" ARM_ERR: R numeric expected",C_result);	
	double C_Emin;
	XL_readNumCell(XL_Emin,C_Emin,			" ARM_ERR: Emin numeric expected",C_result);	
	double C_Lmax;
	XL_readNumCell(XL_Lmax,C_Lmax,			" ARM_ERR: Lmax numeric expected",C_result);	
	double C_gamma0;
	XL_readNumCell(XL_gamma0,C_gamma0,			" ARM_ERR: gamma0 numeric expected",C_result);	
	double C_gamma1;
	XL_readNumCell(XL_gamma1,C_gamma1,			" ARM_ERR: gamma1 numeric expected",C_result);	
	double C_sig;
	XL_readNumCell(XL_sig,C_sig,			" ARM_ERR: sig numeric expected",C_result);	
	double C_lambda;
	XL_readNumCell(XL_lambda,C_lambda,		" ARM_ERR: lambda numeric expected",C_result);	
	double C_sigJ;
	XL_readNumCell(XL_sigJ,C_sigJ,			" ARM_ERR: sigJ numeric expected",C_result);
	double C_r;
	XL_readNumCell(XL_r,C_r,			" ARM_ERR: r numeric expected",C_result);	
	double C_s;
	XL_readNumCell(XL_s,C_s,			" ARM_ERR: s numeric expected",C_result);	
	double C_mu;
	XL_readNumCell(XL_mu,C_mu,			" ARM_ERR: mu numeric expected",C_result);	
	double C_fees;
	XL_readNumCell(XL_fees,C_fees,		" ARM_ERR: fees numeric expected",C_result);	
	double C_volDrift;
	XL_readNumCell(XL_volDrift,C_volDrift,		" ARM_ERR: volDrift numeric expected",C_result);	
	double C_volVol;
	XL_readNumCell(XL_volVol,C_volVol,		" ARM_ERR: volVol numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	vector<double> C_params;
	XL_readNumVector(XL_params,C_params," ARM_ERR: params  Vector: array of numeric expected",C_result);

	/// call the function
	long retCode=ARMLOCAL_JumpDiffusion_Mepi_Call(
		C_P0,
		C_f0,
		C_T,
		C_K,
		C_R,
		C_Emin,
		C_Lmax,
		C_gamma0,
		C_gamma1,
		C_sig,
		C_lambda,
		C_sigJ,
		C_r,
		C_s,
		C_mu,
		C_fees,
		C_volDrift,
		C_volVol,
		C_callput,
		C_params,
		C_result);
	
	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		int nbRows = 2;
		int nbCols = 1;
		VECTOR<double> vectorResult(nbRows*nbCols);
		for ( size_t i=0; i < nbRows; ++i)
			for ( size_t j=0; j < nbCols; ++j)
				vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_JumpDiffusion_Mepi_Option" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_Util_TrigonalSolve(
	
	LPXLOPER XL_A_Vec,
	LPXLOPER XL_B_Vec,
	LPXLOPER XL_C_Vec,
	LPXLOPER XL_R_Vec
	)
{
	ADD_LOG("Local_Util_TrigonalSolve");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_A_Vec;
		XL_readNumVector(XL_A_Vec,C_A_Vec," ARM_ERR: A Vector: array of numeric expected",C_result);
		vector<double> C_B_Vec;
		XL_readNumVector(XL_B_Vec,C_B_Vec," ARM_ERR: B Vector: array of numeric expected",C_result);
		vector<double> C_C_Vec;
		XL_readNumVector(XL_C_Vec,C_C_Vec," ARM_ERR: C Vector: array of numeric expected",C_result);
		vector<double> C_R_Vec;
		XL_readNumVector(XL_R_Vec,C_R_Vec," ARM_ERR: R Vector: array of numeric expected",C_result);

		/// call the function
		long retCode=ARMLOCAL_Util_TrigonalSolve(
			C_A_Vec,
			C_B_Vec,
			C_C_Vec,
			C_R_Vec,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = C_B_Vec.size();
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Util_TrigonalSolve" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR Vanilla Spreadoption Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	


	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_SpreadOption(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_K,
		C_T,
		C_CallPut,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric_Appell(
	LPXLOPER XL_a,
	LPXLOPER XL_b1,
	LPXLOPER XL_b2,
	LPXLOPER XL_c,
	LPXLOPER XL_x,
	LPXLOPER XL_xim,
	LPXLOPER XL_y,
	LPXLOPER XL_yim,
	LPXLOPER XL_nb
	)
		{
			ADD_LOG("Local_Hypergeometric_Appell");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_a;
	XL_readNumCell(XL_a,C_a,				" ARM_ERR: a numeric expected",C_result);
	double C_b1;
	XL_readNumCell(XL_b1,C_b1,				" ARM_ERR: b1 numeric expected",C_result);	
	double C_b2;
	XL_readNumCell(XL_b2,C_b2,				" ARM_ERR: b2 numeric expected",C_result);	
	double C_c;
	XL_readNumCell(XL_c,C_c,				" ARM_ERR: c numeric expected",C_result);	
	double C_x;
	XL_readNumCell(XL_x,C_x,				" ARM_ERR: x numeric expected",C_result);	
	double C_y;
	XL_readNumCell(XL_y,C_y,				" ARM_ERR: y numeric expected",C_result);	
	double C_xim;
	XL_readNumCell(XL_xim,C_xim,				" ARM_ERR: xim numeric expected",C_result);	
	double C_yim;
	XL_readNumCell(XL_yim,C_yim,				" ARM_ERR: yim numeric expected",C_result);	
	double C_nb;
	XL_readNumCell(XL_nb,C_nb,			" ARM_ERR: nb numeric expected",C_result);	
	


	/// call the function
	long retCode=ARMLOCAL_Hypergeometric_Appell(
		C_a,
		C_b1,
		C_b2,
		C_c,
		C_x,
		C_xim,
		C_y,
		C_yim,
		C_nb,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Hypergeometric_Appell" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Util_BiSABR_Eigenvalues(
	
	LPXLOPER XL_rho1,
	LPXLOPER XL_rho2,
	LPXLOPER XL_rhoS,
	LPXLOPER XL_rhoV,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21
	)
{
	ADD_LOG("Local_Util_BiSABR_Eigenvalues");
	return Local_Util_Eigenvalues4(
	XL_rho1,
	XL_rhoS,
	XL_rhoc12,
	XL_rhoc21,
	XL_rhoV,
	XL_rho2);

}

__declspec(dllexport) LPXLOPER WINAPI Local_Util_Eigenvalues4(
	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho14,
	LPXLOPER XL_rho23,
	LPXLOPER XL_rho24,
	LPXLOPER XL_rho34
	)
{
	ADD_LOG("Local_Util_Eigenvalues4");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,				" ARM_ERR: rho12 numeric expected",C_result);	
	double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,				" ARM_ERR: rho13 numeric expected",C_result);	
	double C_rho14;
	XL_readNumCell(XL_rho14,C_rho14,				" ARM_ERR: rho14 numeric expected",C_result);	
	double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,				" ARM_ERR: rho23 numeric expected",C_result);	
	double C_rho24;
	XL_readNumCell(XL_rho24,C_rho24,				" ARM_ERR: rho24 numeric expected",C_result);	
	double C_rho34;
	XL_readNumCell(XL_rho34,C_rho34,				" ARM_ERR: rho34 numeric expected",C_result);	

		/// call the function
		long retCode=ARMLOCAL_Util_Eigenvalues4(
			C_rho12,
			C_rho13,
			C_rho14,
			C_rho23,
			C_rho24,
			C_rho34,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 4;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Util_Eigenvalues4" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Util_Eigenvalues3(
	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23
	)
{
	ADD_LOG("Local_Util_Eigenvalues3");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
	double C_rho12;
	XL_readNumCell(XL_rho12,C_rho12,				" ARM_ERR: rho12 numeric expected",C_result);	
	double C_rho13;
	XL_readNumCell(XL_rho13,C_rho13,				" ARM_ERR: rho13 numeric expected",C_result);	
	double C_rho23;
	XL_readNumCell(XL_rho23,C_rho23,				" ARM_ERR: rho23 numeric expected",C_result);	
		/// call the function
		long retCode=ARMLOCAL_Util_Eigenvalues3(
			C_rho12,
			C_rho13,
			C_rho23,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 3;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Util_Eigenvalues3" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Util_BiSABR_CorrelationEvolution(
	
	LPXLOPER XL_rho1,
	LPXLOPER XL_rho2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_newrho1,
	LPXLOPER XL_newrho2
	)
{
	ADD_LOG("Local_Util_BiSABR_CorrelationEvolution");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,				" ARM_ERR: rho2 numeric expected",C_result);	
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,				" ARM_ERR: rhos numeric expected",C_result);	
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,				" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,				" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,				" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_newrho1;
	XL_readNumCell(XL_newrho1,C_newrho1,				" ARM_ERR: newrho1 numeric expected",C_result);	
	double C_newrho2;
	XL_readNumCell(XL_newrho2,C_newrho2,				" ARM_ERR: newrho2 numeric expected",C_result);	

		/// call the function
		long retCode=ARMLOCAL_Util_BiSABR_CorrelationEvolution(
			C_rho1,
			C_rho2,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_newrho1,
			C_newrho2,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 3;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Util_BiSABR_CorrelationEvolution" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Calibrate(
	
	LPXLOPER	XL_F1_vec, 
	LPXLOPER	XL_alpha1_vec,
	LPXLOPER	XL_beta1_vec ,
	LPXLOPER	XL_rho1_vec,
	LPXLOPER	XL_nu1_vec,
	LPXLOPER	XL_F2_vec,
	LPXLOPER	XL_alpha2_vec,
	LPXLOPER	XL_beta2_vec ,
	LPXLOPER	XL_rho2_vec,
	LPXLOPER	XL_nu2_vec,
	LPXLOPER	XL_strike_vec,
	LPXLOPER	XL_maturity_vec,
	LPXLOPER	XL_callput_vec,
	LPXLOPER	XL_price_vec,
	LPXLOPER	XL_weight_vec,
	LPXLOPER	XL_initialparams

	)

{
	ADD_LOG("Local_BiSABR_Calibrate");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
		vector<double> C_F1_vec;
		XL_readNumVector(XL_F1_vec,C_F1_vec," ARM_ERR: F1_vec Vector: array of numeric expected",C_result);
		vector<double> C_alpha1_vec;
		XL_readNumVector(XL_alpha1_vec,C_alpha1_vec," ARM_ERR: alpha1_vec Vector: array of numeric expected",C_result);
		vector<double> C_beta1_vec;
		XL_readNumVector(XL_beta1_vec,C_beta1_vec," ARM_ERR: beta1_vec Vector: array of numeric expected",C_result);
		vector<double> C_rho1_vec;
		XL_readNumVector(XL_rho1_vec,C_rho1_vec," ARM_ERR: rho1_vec Vector: array of numeric expected",C_result);
		vector<double> C_nu1_vec;
		XL_readNumVector(XL_nu1_vec,C_nu1_vec," ARM_ERR: nu1_vec Vector: array of numeric expected",C_result);

		vector<double> C_F2_vec;
		XL_readNumVector(XL_F2_vec,C_F2_vec," ARM_ERR: F2_vec Vector: array of numeric expected",C_result);
		vector<double> C_alpha2_vec;
		XL_readNumVector(XL_alpha2_vec,C_alpha2_vec," ARM_ERR: alpha2_vec Vector: array of numeric expected",C_result);
		vector<double> C_beta2_vec;
		XL_readNumVector(XL_beta2_vec,C_beta2_vec," ARM_ERR: beta2_vec Vector: array of numeric expected",C_result);
		vector<double> C_rho2_vec;
		XL_readNumVector(XL_rho2_vec,C_rho2_vec," ARM_ERR: rho2_vec Vector: array of numeric expected",C_result);
		vector<double> C_nu2_vec;
		XL_readNumVector(XL_nu2_vec,C_nu2_vec," ARM_ERR: nu2_vec Vector: array of numeric expected",C_result);

		vector<double> C_strike_vec;
		XL_readNumVector(XL_strike_vec,C_strike_vec," ARM_ERR: strike_vec Vector: array of numeric expected",C_result);
		vector<double> C_maturity_vec;
		XL_readNumVector(XL_maturity_vec,C_maturity_vec," ARM_ERR: maturity_vec Vector: array of numeric expected",C_result);
		vector<double> C_callput_vec;
		XL_readNumVector(XL_callput_vec,C_callput_vec," ARM_ERR: callput_vec Vector: array of numeric expected",C_result);
		vector<double> C_price_vec;
		XL_readNumVector(XL_price_vec,C_price_vec," ARM_ERR: price_vec Vector: array of numeric expected",C_result);

		vector<double> C_weight_vec;
		XL_readNumVector(XL_weight_vec,C_weight_vec," ARM_ERR: weight_vec Vector: array of numeric expected",C_result);

	
		vector<double> C_initialparams;
		XL_readNumVector(XL_initialparams,C_initialparams," ARM_ERR: initialparams Vector: array of numeric expected",C_result);




		/// call the function
		long retCode=ARMLOCAL_BiSABR_Calibrate(
			C_F1_vec, 
			C_alpha1_vec,
			C_beta1_vec ,
			C_rho1_vec,
			C_nu1_vec,
			C_F2_vec,
			C_alpha2_vec,
			C_beta2_vec ,
			C_rho2_vec,
			C_nu2_vec,
			C_strike_vec,
			C_maturity_vec,
			C_callput_vec,
			C_price_vec,
			C_weight_vec,
			C_initialparams,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = C_F1_vec.size();
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///	Lognormal Vanilla digital option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_LN_DigitalOption(
	LPXLOPER XL_forward,
	LPXLOPER XL_strike,
	LPXLOPER XL_maturity,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_volatility
	)
		{
			ADD_LOG("Local_LN_DigitalOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_forward;
	XL_readNumCell(XL_forward,C_forward,				" ARM_ERR: forward numeric expected",C_result);
	double C_strike;
	XL_readNumCell(XL_strike,C_strike,				" ARM_ERR: strike numeric expected",C_result);	
	double C_maturity;
	XL_readNumCell(XL_maturity,C_maturity,				" ARM_ERR: maturity numeric expected",C_result);	
	double C_volatility;
	XL_readNumCell(XL_volatility,C_volatility,				" ARM_ERR: volatility numeric expected",C_result);	
		
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	


	/// call the function
	long retCode=ARMLOCAL_LN_DigitalOption(
		C_forward,
		C_strike,
		C_maturity,
		C_CallPut,
		C_volatility,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LN_DigitalOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	Lognormal Ratio option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_LN_RatioOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_Mu1,
	LPXLOPER XL_Sigma1,
	LPXLOPER XL_S2,
	LPXLOPER XL_Mu2,
	LPXLOPER XL_Sigma2,
	LPXLOPER XL_Rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut
	)
		{
			ADD_LOG("Local_LN_RatioOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,				" ARM_ERR: S1 numeric expected",C_result);
	double C_Mu1;
	XL_readNumCell(XL_Mu1,C_Mu1,				" ARM_ERR: Mu1 numeric expected",C_result);	
	double C_Sigma1;
	XL_readNumCell(XL_Sigma1,C_Sigma1,				" ARM_ERR: Sigma1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,				" ARM_ERR: S2 numeric expected",C_result);	
	double C_Mu2;
	XL_readNumCell(XL_Mu2,C_Mu2,				" ARM_ERR: Mu2 numeric expected",C_result);
	double C_Sigma2;
	XL_readNumCell(XL_Sigma2,C_Sigma2,				" ARM_ERR: Sigma2 numeric expected",C_result);	
	double C_Rho;
	XL_readNumCell(XL_Rho,C_Rho,				" ARM_ERR: Rho numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);	
	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	


	/// call the function
	long retCode=ARMLOCAL_LN_RatioOption(
		C_S1,
		C_Mu1,
		C_Sigma1,
		C_S2,
		C_Mu2,
		C_Sigma2,
		C_Rho,
		C_K,
		C_T,
		C_CallPut,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LN_RatioOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	Lognormal Product option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_LN_ProductOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_Mu1,
	LPXLOPER XL_Sigma1,
	LPXLOPER XL_S2,
	LPXLOPER XL_Mu2,
	LPXLOPER XL_Sigma2,
	LPXLOPER XL_Rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut
	)
		{
			ADD_LOG("Local_LN_ProductOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_S1;
	XL_readNumCell(XL_S1,C_S1,				" ARM_ERR: S1 numeric expected",C_result);
	double C_Mu1;
	XL_readNumCell(XL_Mu1,C_Mu1,				" ARM_ERR: Mu1 numeric expected",C_result);	
	double C_Sigma1;
	XL_readNumCell(XL_Sigma1,C_Sigma1,				" ARM_ERR: Sigma1 numeric expected",C_result);	
	double C_S2;
	XL_readNumCell(XL_S2,C_S2,				" ARM_ERR: S2 numeric expected",C_result);	
	double C_Mu2;
	XL_readNumCell(XL_Mu2,C_Mu2,				" ARM_ERR: Mu2 numeric expected",C_result);
	double C_Sigma2;
	XL_readNumCell(XL_Sigma2,C_Sigma2,				" ARM_ERR: Sigma2 numeric expected",C_result);	
	double C_Rho;
	XL_readNumCell(XL_Rho,C_Rho,				" ARM_ERR: Rho numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);	
	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	


	/// call the function
	long retCode=ARMLOCAL_LN_ProductOption(
		C_S1,
		C_Mu1,
		C_Sigma1,
		C_S2,
		C_Mu2,
		C_Sigma2,
		C_Rho,
		C_K,
		C_T,
		C_CallPut,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LN_ProductOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	copula base SABR calls
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCall(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_LegendreNb,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh
	)


		{
			ADD_LOG("Local_SABR_GaussianSABRDigitalCall");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,				" ARM_ERR: f1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);
	double C_flag1;
	XL_GETCONVSABRFLAG(XL_flag1,C_flag1," ARM_ERR: SABR flag1 string expected",C_result);
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,				" ARM_ERR: f2 numeric expected",C_result);
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,				" ARM_ERR: alpha2 numeric expected",C_result);
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,				" ARM_ERR: beta2 numeric expected",C_result);
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,				" ARM_ERR: rho2 numeric expected",C_result);
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,				" ARM_ERR: nu2 numeric expected",C_result);
	double C_flag2;
	XL_GETCONVSABRFLAG(XL_flag2,C_flag2," ARM_ERR: SABR flag2 string expected",C_result);
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,				" ARM_ERR: rho numeric expected",C_result);

	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);

	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);

	double C_LegendreNb;
	double LegendreNb_default=120;
    XL_readNumCellWD(XL_LegendreNb, C_LegendreNb, LegendreNb_default, " ARM_ERR: LegendreNb: numeric expected",C_result);
	double C_alpha_exp;
	double alpha_exp_default=0.01;
    XL_readNumCellWD(XL_alpha_exp, C_alpha_exp, alpha_exp_default, " ARM_ERR: alpha_exp: numeric expected",C_result);
	double C_alpha_tanh;
	double alpha_tanh_default=1.5;
    XL_readNumCellWD(XL_alpha_tanh, C_alpha_tanh, alpha_tanh_default, " ARM_ERR: alpha_tanh: numeric expected",C_result);
	double C_kb_tanh;
	double kb_tanh_default=0.02;
    XL_readNumCellWD(XL_kb_tanh, C_kb_tanh, kb_tanh_default, " ARM_ERR: kb_tanh: numeric expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_SABR_GaussianSABRDigitalCall(
		C_f1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_flag1,
		C_f2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_flag2,
		C_rho,
		C_K,
		C_T,
		C_LegendreNb,C_alpha_exp,C_alpha_tanh,C_kb_tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_GaussianSABRDigitalCall" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCallPayingS1(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_LegendreNb,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh
	)


		{
			ADD_LOG("Local_SABR_GaussianSABRDigitalCallPayingS1");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,				" ARM_ERR: f1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);
	double C_flag1;
	XL_GETCONVSABRFLAG(XL_flag1,C_flag1," ARM_ERR: SABR flag1 string expected",C_result);
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,				" ARM_ERR: f2 numeric expected",C_result);
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,				" ARM_ERR: alpha2 numeric expected",C_result);
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,				" ARM_ERR: beta2 numeric expected",C_result);
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,				" ARM_ERR: rho2 numeric expected",C_result);
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,				" ARM_ERR: nu2 numeric expected",C_result);
	double C_flag2;
	XL_GETCONVSABRFLAG(XL_flag2,C_flag2," ARM_ERR: SABR flag2 string expected",C_result);
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,				" ARM_ERR: rho numeric expected",C_result);

	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);

	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);

	double C_LegendreNb;
	double LegendreNb_default=120;
    XL_readNumCellWD(XL_LegendreNb, C_LegendreNb, LegendreNb_default, " ARM_ERR: LegendreNb: numeric expected",C_result);
	double C_alpha_exp;
	double alpha_exp_default=0.01;
    XL_readNumCellWD(XL_alpha_exp, C_alpha_exp, alpha_exp_default, " ARM_ERR: alpha_exp: numeric expected",C_result);
	double C_alpha_tanh;
	double alpha_tanh_default=1.5;
    XL_readNumCellWD(XL_alpha_tanh, C_alpha_tanh, alpha_tanh_default, " ARM_ERR: alpha_tanh: numeric expected",C_result);
	double C_kb_tanh;
	double kb_tanh_default=0.02;
    XL_readNumCellWD(XL_kb_tanh, C_kb_tanh, kb_tanh_default, " ARM_ERR: kb_tanh: numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS1(
		C_f1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_flag1,
		C_f2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_flag2,
		C_rho,
		C_K,
		C_T,
		C_LegendreNb,C_alpha_exp,C_alpha_tanh,C_kb_tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_GaussianSABRDigitalCallPayingS1" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCallPayingS2(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_LegendreNb,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh
	)


		{
			ADD_LOG("Local_SABR_GaussianSABRDigitalCallPayingS2");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,				" ARM_ERR: f1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);
	double C_flag1;
	XL_GETCONVSABRFLAG(XL_flag1,C_flag1," ARM_ERR: SABR flag1 string expected",C_result);
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,				" ARM_ERR: f2 numeric expected",C_result);
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,				" ARM_ERR: alpha2 numeric expected",C_result);
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,				" ARM_ERR: beta2 numeric expected",C_result);
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,				" ARM_ERR: rho2 numeric expected",C_result);
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,				" ARM_ERR: nu2 numeric expected",C_result);
	double C_flag2;
	XL_GETCONVSABRFLAG(XL_flag2,C_flag2," ARM_ERR: SABR flag2 string expected",C_result);
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,				" ARM_ERR: rho numeric expected",C_result);

	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);

	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);

	double C_LegendreNb;
	double LegendreNb_default=120;
    XL_readNumCellWD(XL_LegendreNb, C_LegendreNb, LegendreNb_default, " ARM_ERR: LegendreNb: numeric expected",C_result);
	double C_alpha_exp;
	double alpha_exp_default=0.01;
    XL_readNumCellWD(XL_alpha_exp, C_alpha_exp, alpha_exp_default, " ARM_ERR: alpha_exp: numeric expected",C_result);
	double C_alpha_tanh;
	double alpha_tanh_default=1.5;
    XL_readNumCellWD(XL_alpha_tanh, C_alpha_tanh, alpha_tanh_default, " ARM_ERR: alpha_tanh: numeric expected",C_result);
	double C_kb_tanh;
	double kb_tanh_default=0.02;
    XL_readNumCellWD(XL_kb_tanh, C_kb_tanh, kb_tanh_default, " ARM_ERR: kb_tanh: numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS2(
		C_f1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_flag1,
		C_f2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_flag2,
		C_rho,
		C_K,
		C_T,
		C_LegendreNb,C_alpha_exp,C_alpha_tanh,C_kb_tanh,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_GaussianSABRDigitalCallPayingS2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCallPayingS3(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_f3,
	LPXLOPER XL_sigma3,
	LPXLOPER XL_Correlations_Vec,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Params_Vec

	)


		{
			ADD_LOG("Local_SABR_GaussianSABRDigitalCallPayingS3");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,				" ARM_ERR: f1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);
	double C_flag1;
	XL_GETCONVSABRFLAG(XL_flag1,C_flag1," ARM_ERR: SABR flag1 string expected",C_result);
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,				" ARM_ERR: f2 numeric expected",C_result);
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,				" ARM_ERR: alpha2 numeric expected",C_result);
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,				" ARM_ERR: beta2 numeric expected",C_result);
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,				" ARM_ERR: rho2 numeric expected",C_result);
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,				" ARM_ERR: nu2 numeric expected",C_result);
	double C_flag2;
	XL_GETCONVSABRFLAG(XL_flag2,C_flag2," ARM_ERR: SABR flag2 string expected",C_result);
	double C_f3;
	XL_readNumCell(XL_f3,C_f3,				" ARM_ERR: f3 numeric expected",C_result);
	double C_sigma3;
	XL_readNumCell(XL_sigma3,C_sigma3,				" ARM_ERR: sigma3 numeric expected",C_result);

	VECTOR<double> C_Correlations_Vec;
	XL_readNumVector(XL_Correlations_Vec,C_Correlations_Vec," ARM_ERR: Correlations_Vec: array of numeric expected",C_result);


	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);
	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);


	VECTOR<double> C_Params_Vec;
	XL_readNumVector(XL_Params_Vec,C_Params_Vec," ARM_ERR: Params_Vec: array of numeric expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS3(
		C_f1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_flag1,
		C_f2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_flag2,
		C_f3,
		C_sigma3,
		C_Correlations_Vec,
		C_K,
		C_T,
		C_Params_Vec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS3" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_TarnProxyCommon(
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Fwds,
	LPXLOPER XL_DensityFunctors,
	LPXLOPER XL_DiscountFactors,
	LPXLOPER XL_LevPrec,	
	LPXLOPER XL_Lev,		
	LPXLOPER XL_Fix,
	LPXLOPER XL_Cap,
	LPXLOPER XL_Floor,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Dcf,
	LPXLOPER XL_Target,
	LPXLOPER XL_GlobalCap,
	LPXLOPER XL_GlobalFloor,
	LPXLOPER XL_CorrelInput,
	LPXLOPER XL_NbSimul,
	bool PersistentInXL )
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		vector<double> defaultVector(0);
		defaultVector.clear();

		vector<double> C_ResetDates;
		XL_readNumVectorWD(XL_ResetDates,C_ResetDates,defaultVector," ARM_ERR: ResetDates: array of numeric expected",C_result);

		vector<double> C_Fwds;
		vector<double> C_DiscountFactors;
		vector<double> C_LevPrec;
		vector<double> C_Lev;	
		vector<double> C_Fix;
		vector<double> C_Cap;
		vector<double> C_Floor;
		vector<double> C_Fees;
		vector<double> C_Dcf;
		if(C_ResetDates.size() > 0)
        {
			XL_readNumVector(XL_Fwds,C_Fwds," ARM_ERR: Fwds: array of numeric expected",C_result);
			XL_readNumVector(XL_DiscountFactors,C_DiscountFactors," ARM_ERR: DiscountFactors: array of numeric expected",C_result);
			XL_readNumVector(XL_LevPrec,C_LevPrec," ARM_ERR: LevPrec: array of numeric expected",C_result);
			XL_readNumVector(XL_Lev,C_Lev," ARM_ERR: Lev: array of numeric expected",C_result);
			XL_readNumVector(XL_Fix,C_Fix," ARM_ERR: Fix: array of numeric expected",C_result);
			XL_readNumVector(XL_Cap,C_Cap," ARM_ERR: Cap: array of numeric expected",C_result);
			XL_readNumVector(XL_Floor,C_Floor," ARM_ERR: Floor: array of numeric expected",C_result);
			XL_readNumVector(XL_Fees,C_Fees," ARM_ERR: Fees: array of numeric expected",C_result);
			XL_readNumVector(XL_Dcf,C_Dcf," ARM_ERR: Fees: array of numeric expected",C_result);
        }

		VECTOR<CCString> C_DensityFunctorId;
		XL_readStrVector (XL_DensityFunctors,C_DensityFunctorId," ARM_ERR: DensityFunctorId: array of object expected",DOUBLE_TYPE,C_result);

		size_t size = C_DensityFunctorId.size();
		VECTOR<long> C_DensityFunctorIdLong (size);
		size_t i;
		for(i = 0; i < size; ++i )
		{
			C_DensityFunctorIdLong[i] = LocalGetNumObjectId(C_DensityFunctorId[i]); 
		}

		double C_Target;
		XL_readNumCell( XL_Target, C_Target, " ARM_ERR: Target: numeric expected",C_result);
		
		CCString C_GlobalCap;
		XL_readStrCellWD( XL_GlobalCap,C_GlobalCap,"Y"," ARM_ERR: GlobalCap: string expected",C_result);
		bool C_GlobalCapBool;
		C_GlobalCap.toUpper();
		if (C_GlobalCap == "Y" )
			C_GlobalCapBool = true;
		else if (C_GlobalCap == "N")
			C_GlobalCapBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for GlobalCap");
		
		CCString C_GlobalFloor;
		XL_readStrCellWD( XL_GlobalFloor,C_GlobalFloor,"Y"," ARM_ERR: GlobalFloor: string expected",C_result);
		bool C_GlobalFloorBool;
		C_GlobalFloor.toUpper();
		if (C_GlobalFloor == "Y" )
			C_GlobalFloorBool = true;
		else if (C_GlobalFloor == "N")
			C_GlobalFloorBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for GlobalFloor");

		double C_CorrelInput;
		XL_readNumCell( XL_CorrelInput, C_CorrelInput, " ARM_ERR: CorrelInput: numeric expected",C_result);
		
		double C_NbSimul;
		XL_readNumCell( XL_NbSimul, C_NbSimul, " ARM_ERR: NbSimul: numeric expected",C_result);
		int C_NbSimulInt = C_NbSimul;

		exportFunc16Args<	vector<double>, 
							vector<double>,
							VECTOR<long>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							double,
							bool,
							bool,
							double,
							int
						> 
			ourFunc(
				C_ResetDates,
				C_Fwds,
				C_DensityFunctorIdLong,
				C_DiscountFactors,
				C_LevPrec,
				C_Lev,
				C_Fix,
				C_Cap,
				C_Floor,
				C_Fees,
				C_Dcf,
				C_Target,
				C_GlobalCapBool,
				C_GlobalFloorBool,
				C_CorrelInput,
				C_NbSimulInt,
				ARMLOCAL_TarnProxy_Create);

		/// call the general function
		fillXL_Result( LOCAL_PROXY_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TarnProxy_Create(
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Fwds,
	LPXLOPER XL_DensityFunctors,
	LPXLOPER XL_DiscountFactors,
	LPXLOPER XL_LevPrec,	
	LPXLOPER XL_Lev,		
	LPXLOPER XL_Fix,
	LPXLOPER XL_Cap,
	LPXLOPER XL_Floor,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Dcf,
	LPXLOPER XL_Target,
	LPXLOPER XL_GlobalCap,
	LPXLOPER XL_GlobalFloor,
	LPXLOPER XL_CorrelInput,
	LPXLOPER XL_NbSimul)
{
	ADD_LOG("Local_TarnProxy_Create");
	bool PersistentInXL = true;
	return Local_TarnProxyCommon(
		XL_ResetDates,
		XL_Fwds,
		XL_DensityFunctors,
		XL_DiscountFactors,
		XL_LevPrec,	
		XL_Lev,		
		XL_Fix,
		XL_Cap,
		XL_Floor,
		XL_Fees,
		XL_Dcf,
		XL_Target,
		XL_GlobalCap,
		XL_GlobalFloor,
		XL_CorrelInput,
		XL_NbSimul,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TarnProxy_Create(
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Fwds,
	LPXLOPER XL_DensityFunctors,
	LPXLOPER XL_DiscountFactors,
	LPXLOPER XL_LevPrec,	
	LPXLOPER XL_Lev,		
	LPXLOPER XL_Fix,
	LPXLOPER XL_Cap,
	LPXLOPER XL_Floor,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Dcf,
	LPXLOPER XL_Target,
	LPXLOPER XL_GlobalCap,
	LPXLOPER XL_GlobalFloor,
	LPXLOPER XL_CorrelInput,
	LPXLOPER XL_NbSimul)
{
	ADD_LOG("Local_PXL_TarnProxy_Create");
	bool PersistentInXL = false;
	return Local_TarnProxyCommon(
		XL_ResetDates,
		XL_Fwds,
		XL_DensityFunctors,
		XL_DiscountFactors,
		XL_LevPrec,	
		XL_Lev,		
		XL_Fix,
		XL_Cap,
		XL_Floor,
		XL_Fees,
		XL_Dcf,
		XL_Target,
		XL_GlobalCap,
		XL_GlobalFloor,
		XL_CorrelInput,
		XL_NbSimul,
		PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_TarnProxy_GetPrice(
	LPXLOPER XL_TarnProxyId,
    LPXLOPER XL_NbPayoff)
{
	ADD_LOG("Local_TarnProxy_GetPrice");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		double C_NbPayoffDef=0;
		double C_NbPayoff;
		XL_readNumCellWD( XL_NbPayoff, C_NbPayoff,C_NbPayoffDef, " ARM_ERR: NbPayoff: numeric expected",C_result);
		int C_NbPayoffInt = C_NbPayoff;

		CCString C_TarnProxyStrId;
		XL_readStrCell( XL_TarnProxyId, C_TarnProxyStrId,	" ARM_ERR: TarnProxy Id: Object expected", C_result);
		long C_TarnProxyId = LocalGetNumObjectId(C_TarnProxyStrId);


		long retCode = ARMLOCAL_TarnProxy_GetPrice(
				C_TarnProxyId,
				C_NbPayoffInt,
				C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_VBMinMaxProxy_Common(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_totalVol,
	LPXLOPER XL_leftVol,
	LPXLOPER XL_rightVol,
	LPXLOPER XL_nu,
	LPXLOPER XL_rho,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_sabrDiff,
	LPXLOPER XL_typeprice,
	LPXLOPER XL_type1sens,
	LPXLOPER XL_type2sens,
	LPXLOPER XL_rate1Lev,
	LPXLOPER XL_rate1Add,
	LPXLOPER XL_capRateLev,
	LPXLOPER XL_capRateAdd,
	LPXLOPER XL_vbLev,
	LPXLOPER XL_maxChoice,
	LPXLOPER XL_minChoice,
	LPXLOPER XL_minmaxFreq,
	bool PersistentInXL)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		vector<double> defaultVector(0);
		defaultVector.clear();

		double C_AsOf;
		vector<double>	C_resetDates, C_fwdRates, C_totalVol, C_leftVol, C_rightVol, C_nu, C_rho;
		double C_nbSimul;
		double C_typeprice;
		int C_type1sens, C_type2sens, C_maxChoice, C_minChoice;
		bool C_sabrDiff;
		CCString sabrdiff, sens1, sens2, maxchoice, minchoice;
		double C_minmaxFreq;
		vector<double> C_rate1Lev, C_rate1Add, C_capRateLev, C_capRateAdd, C_vbLev;

		XL_readNumCell(XL_AsOf, C_AsOf, "ARM_ERR: asof numeric expected",C_result);

		XL_readNumVector(XL_resetDates, C_resetDates, "ARM_ERR : ResetDates : array of numeric expected", C_result);
		XL_readNumVector(XL_fwdRates, C_fwdRates, "ARM_ERR : fwdRates : array of numeric expected", C_result);
		XL_readNumVector(XL_totalVol, C_totalVol, "ARM_ERR : Spot Vol : array of numeric expected", C_result);
		XL_readNumVector(XL_leftVol, C_leftVol, "ARM_ERR : left Vol : array of numeric expected", C_result);
		XL_readNumVector(XL_rightVol, C_rightVol, "ARM_ERR : right Vol : array of numeric expected", C_result);
		XL_readNumVector(XL_nu, C_nu, "ARM_ERR : nu : array of numeric expected", C_result);
		XL_readNumVector(XL_rho, C_rho, "ARM_ERR : rho : array of numeric expected", C_result);

		XL_readNumCell(XL_nbSimul, C_nbSimul, "ARM_ERR : nb simul : numeric expected", C_result);

		XL_readStrCellWD(XL_sabrDiff,sabrdiff,"Y","ARM_ERR : sabr diff : string expected", C_result);
		sabrdiff.toUpper();
		C_sabrDiff = sabrdiff == "Y" ? true : false;

		XL_readNumCell(XL_typeprice, C_typeprice," ARM_ERR : type : numeric expected",C_result);

		XL_readStrCellWD(XL_type1sens, sens1, "C"," ARM_ERR : sens 1 : string expected",C_result);
		XL_readStrCellWD(XL_type2sens, sens2, "C"," ARM_ERR : sens 2 : string expected",C_result);

		sens1.toUpper();
		sens2.toUpper();

		C_type1sens = sens1 == "CALL" || sens1 == "C" ? 1 : sens1 == "PUT" || sens1 == "P" ? -1 : 0;
		C_type2sens = sens2 == "CALL" || sens2 == "C" ? 1 : sens2 == "PUT" || sens2 == "P" ? -1 : 0;

		XL_readNumVectorWD(XL_rate1Lev, C_rate1Lev, defaultVector, "ARM_ERR : first rate leverage : array of numeric expected", C_result);
		XL_readNumVectorWD(XL_rate1Add, C_rate1Add, defaultVector, "ARM_ERR : first rate add : array of numeric expected", C_result);
		XL_readNumVectorWD(XL_capRateLev, C_capRateLev, defaultVector, "ARM_ERR : cap first rate leverage : array of numeric expected", C_result);
		XL_readNumVectorWD(XL_capRateAdd, C_capRateAdd, defaultVector, "ARM_ERR : cap first rate add : array of numeric expected", C_result);
		XL_readNumVectorWD(XL_vbLev, C_vbLev, defaultVector, "ARM_ERR : std vol bond leverage : array of numeric expected", C_result);

		if(C_rate1Lev.size() == 0) C_rate1Lev.resize(C_resetDates.size(), 1.);
		if(C_rate1Add.size() == 0) C_rate1Add.resize(C_resetDates.size(), 0.);
		if(C_capRateLev.size() == 0) C_capRateLev.resize(C_resetDates.size(), 0.);
		if(C_capRateAdd.size() == 0) C_capRateAdd.resize(C_resetDates.size(), 0.);
		if(C_vbLev.size() == 0) C_vbLev.resize(C_resetDates.size(), 0.);

		XL_readStrCellWD( XL_maxChoice, maxchoice,"max", "ARM_ERR : max choice : string expected", C_result);
		XL_readStrCellWD( XL_minChoice, minchoice,"min", "ARM_ERR : min choice : string expected", C_result);
		maxchoice.toUpper();
		minchoice.toUpper();

		C_maxChoice = maxchoice == "DEBUT" ? 0 : maxchoice == "FIN" ? 1 : 2;
		C_minChoice = minchoice == "DEBUT" ? 0 : minchoice == "FIN" ? 1 : 2;

		long resetfreq;
		XL_GETFREQUENCYWD(XL_minmaxFreq, resetfreq, "Z", " ARM_ERR: reset freq: string expected",	C_result );

		switch(resetfreq)
		{
		case K_ANNUAL:
			C_minmaxFreq = 1.;
			break;

		case K_SEMIANNUAL:
			C_minmaxFreq = 0.5;
			break;

		case K_QUARTERLY:
			C_minmaxFreq = 0.25;
			break;

		case K_BIMONTHLY:
			C_minmaxFreq = 1./ 6.;
			break;

		case K_MONTHLY:
			C_minmaxFreq = 1./12.;
			break;

		case K_WEEKLY:
			C_minmaxFreq = 1./52.;
			break;

		case K_DAILY:
			C_minmaxFreq = 1./250.;
			break;

		case K_ZEROCOUPON:
			C_minmaxFreq = 0.;
			break;
		}

		exportFunc21Args<	double,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							int,
							bool,
							int,
							int,
							int,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							int,
							int,
							double
						> 
			ourFunc(
					C_AsOf,
					C_resetDates,
					C_fwdRates,
					C_totalVol,
					C_leftVol,
					C_rightVol,
					C_nu,
					C_rho,
					(int)C_nbSimul,
					C_sabrDiff,
					(int)C_typeprice,
					C_type1sens,
					C_type2sens,
					C_rate1Lev,
					C_rate1Add,
					C_capRateLev,
					C_capRateAdd,
					C_vbLev,
					C_maxChoice,
					C_minChoice,
					C_minmaxFreq,
					ARMLOCAL_VBMinMaxProxy_Create);

		/// call the general function
		fillXL_Result( LOCAL_PROXY_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
	
_declspec(dllexport) LPXLOPER WINAPI Local_PXL_VBMinMaxProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_totalVol,
	LPXLOPER XL_leftVol,
	LPXLOPER XL_rightVol,
	LPXLOPER XL_nu,
	LPXLOPER XL_rho,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_sabrDiff,
	LPXLOPER XL_typeprice,
	LPXLOPER XL_type1sens,
	LPXLOPER XL_type2sens,
	LPXLOPER XL_rate1Lev,
	LPXLOPER XL_rate1Add,
	LPXLOPER XL_capRateLev,
	LPXLOPER XL_capRateAdd,
	LPXLOPER XL_vbLev,
	LPXLOPER XL_maxChoice,
	LPXLOPER XL_minChoice,
	LPXLOPER XL_minmaxFreq
	)
{
	bool PersistentInXL = false;

	return Local_VBMinMaxProxy_Common(XL_AsOf,XL_resetDates, XL_fwdRates, XL_totalVol, XL_leftVol, XL_rightVol,
		XL_nu, XL_rho, XL_nbSimul, XL_sabrDiff, XL_typeprice, XL_type1sens, XL_type2sens, XL_rate1Lev, XL_rate1Add, 
		XL_capRateLev, XL_capRateAdd, XL_vbLev, XL_maxChoice, XL_minChoice, XL_minmaxFreq, PersistentInXL);
}

_declspec(dllexport) LPXLOPER WINAPI Local_VBMinMaxProxy_GetInfo(
	LPXLOPER XL_ProxyId,
	LPXLOPER XL_Info
	)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		CCString C_ProxyStrId;
		XL_readStrCell( XL_ProxyId, C_ProxyStrId,	" ARM_ERR: Proxy Id: Object expected", C_result);
		long C_ProxyId = LocalGetNumObjectId(C_ProxyStrId);
		
		CCString C_InfoStr;
		XL_readStrCell( XL_Info, C_InfoStr, "ARM_ERR : Info : string expected", C_result);
		
		int info;
		C_InfoStr.toUpper();

		if(C_InfoStr == "VB1")
			info = 1;
		else if(C_InfoStr == "VB2")
			info = 2;
		else if(C_InfoStr == "VB3")
			info = 3;
		else if(C_InfoStr == "MAXRATE")
			info = 31;
		else if(C_InfoStr == "MINRATE")
			info = 32;
		else
			info = 0;

		vector<double> vresult;

		long retCode = ARMLOCAL_VBMinMaxProxy_GetInfo(C_ProxyId,info,vresult,C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
			XL_writeNumVector(XL_result, vresult, " ARM_ERR: Could not get result data", C_result);
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_Berm2DatesProxy_GetPrice(
	LPXLOPER XL_ProxyId,
	LPXLOPER XL_Info
	)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		CCString C_ProxyStrId;
		XL_readStrCell( XL_ProxyId, C_ProxyStrId,	" ARM_ERR: Proxy Id: Object expected", C_result);
		long C_ProxyId = LocalGetNumObjectId(C_ProxyStrId);
		
		CCString C_InfoStr;
		XL_readStrCellWD( XL_Info, C_InfoStr, "PRICE", "ARM_ERR : Info : string expected", C_result);
		
		int info;
		C_InfoStr.toUpper();

		if(C_InfoStr == "EUR1")
			info = 1;
		else if(C_InfoStr == "EUR2")
			info = 2;
		else 
			info = 0;

		long retCode = ARMLOCAL_Berm2DatesProxy_GetPrice(C_ProxyId,info,C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


	
_declspec(dllexport) LPXLOPER WINAPI Local_VBMinMaxProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_totalVol,
	LPXLOPER XL_leftVol,
	LPXLOPER XL_rightVol,
	LPXLOPER XL_nu,
	LPXLOPER XL_rho,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_sabrDiff,
	LPXLOPER XL_typeprice,
	LPXLOPER XL_type1sens,
	LPXLOPER XL_type2sens,
	LPXLOPER XL_rate1Lev,
	LPXLOPER XL_rate1Add,
	LPXLOPER XL_capRateLev,
	LPXLOPER XL_capRateAdd,
	LPXLOPER XL_vbLev,
	LPXLOPER XL_maxChoice,
	LPXLOPER XL_minChoice,
	LPXLOPER XL_minmaxFreq
	)
{
	bool PersistentInXL = true;

	return Local_VBMinMaxProxy_Common(XL_AsOf,XL_resetDates, XL_fwdRates, XL_totalVol, XL_leftVol, XL_rightVol,
		XL_nu, XL_rho, XL_nbSimul, XL_sabrDiff, XL_typeprice, XL_type1sens, XL_type2sens, XL_rate1Lev, XL_rate1Add, 
		XL_capRateLev, XL_capRateAdd, XL_vbLev, XL_maxChoice, XL_minChoice, XL_minmaxFreq, PersistentInXL);
}




_declspec(dllexport) LPXLOPER WINAPI Local_SpreadVBProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Fwds1,
	LPXLOPER XL_FwdsVols1,
	LPXLOPER XL_FwdsPartialVols1,
	LPXLOPER XL_Vols1,
	LPXLOPER XL_VVols1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Fwds2,
	LPXLOPER XL_FwdsVols2,
	LPXLOPER XL_FwdsPartialVols2,
	LPXLOPER XL_Vols2,
	LPXLOPER XL_VVols2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_R12Correl,
	LPXLOPER XL_Type,
	LPXLOPER XL_NbSimul,
	LPXLOPER XL_fwdLev,
	LPXLOPER XL_fwdStrikes,
	LPXLOPER XL_Type1Sens,
	LPXLOPER XL_Type2Sens,
	LPXLOPER XL_SabrDiff,
	LPXLOPER XL_CorrMeanRev,
	LPXLOPER XL_CorrVol,
	LPXLOPER XL_FixedCpn,
	LPXLOPER XL_Levier)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		vector<double> defaultVector(0);
		defaultVector.clear();

		vector<double> C_ResetDates;
		XL_readNumVector(XL_ResetDates,C_ResetDates," ARM_ERR: ResetDates: array of numeric expected",C_result);

		vector<double> C_Fwds1;
		XL_readNumVector(XL_Fwds1,C_Fwds1," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_FwdsVols1;
		XL_readNumVector(XL_FwdsVols1,C_FwdsVols1," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_FwdsPartVols1;
		XL_readNumVector(XL_FwdsPartialVols1,C_FwdsPartVols1," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_Vols1;
		XL_readNumVector(XL_Vols1,C_Vols1," ARM_ERR: Vols: array of numeric expected",C_result);
		vector<double> C_VVols1;
		XL_readNumVector(XL_VVols1,C_VVols1," ARM_ERR: Vols: array of numeric expected",C_result);
		vector<double> C_Rho1;
		XL_readNumVector(XL_Rho1,C_Rho1," ARM_ERR: Vols: array of numeric expected",C_result);

		vector<double> C_Fwds2;
		XL_readNumVector(XL_Fwds2,C_Fwds2," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_FwdsVols2;
		XL_readNumVector(XL_FwdsVols2,C_FwdsVols2," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_FwdsPartVols2;
		XL_readNumVector(XL_FwdsPartialVols2,C_FwdsPartVols2," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_Vols2;
		XL_readNumVector(XL_Vols2,C_Vols2," ARM_ERR: Vols: array of numeric expected",C_result);
		vector<double> C_VVols2;
		XL_readNumVector(XL_VVols2,C_VVols2," ARM_ERR: Vols: array of numeric expected",C_result);
		vector<double> C_Rho2;
		XL_readNumVector(XL_Rho2,C_Rho2," ARM_ERR: Vols: array of numeric expected",C_result);

		vector<double> C_R12Correl;
		XL_readNumVector(XL_R12Correl,C_R12Correl," ARM_ERR: Vols: array of numeric expected",C_result);

		vector<double> C_fwdLev, C_fwdStrikes, C_Levier, C_Fixed, defvec;
		XL_readNumVectorWD(XL_fwdLev,C_fwdLev,defvec," ARM_ERR: Vols: array of numeric expected",C_result);
		XL_readNumVectorWD(XL_fwdStrikes,C_fwdStrikes,defvec," ARM_ERR: Vols: array of numeric expected",C_result);
		XL_readNumVectorWD(XL_FixedCpn,C_Fixed,defvec," ARM_ERR: Vols: array of numeric expected",C_result);
		XL_readNumVectorWD(XL_Levier,C_Levier,defvec," ARM_ERR: Vols: array of numeric expected",C_result);
		if(C_fwdLev.size() == 0) C_fwdLev.resize(C_ResetDates.size(),1.);
		if(C_fwdStrikes.size() == 0) C_fwdStrikes.resize(C_ResetDates.size(),0.);
		if(C_Levier.size() == 0) C_Levier.resize(C_ResetDates.size(),0.);
		if(C_Fixed.size() == 0) C_Fixed.resize(C_ResetDates.size(), 0.);

		double C_type;
		XL_readNumCell(XL_Type,C_type," ARM_ERR : type : numeric expected",C_result);
		int C_typevb = (int)C_type;
		double defsens = 1.;
		CCString C_Sens1Str, C_Sens2Str;
		int C_type1sens, C_type2sens;
		XL_readStrCellWD(XL_Type1Sens, C_Sens1Str, "S"," ARM_ERR : type : string expected",C_result);
		XL_readStrCellWD(XL_Type2Sens, C_Sens2Str, "S"," ARM_ERR : type : string expected",C_result);
		C_Sens1Str.toUpper();
		C_Sens2Str.toUpper();
		if(C_Sens1Str == "CALL" || C_Sens1Str == "C")
			C_type1sens = 1;
		else if(C_Sens1Str == "PUT" || C_Sens1Str == "P")
			C_type1sens = -1;
		else
			C_type1sens = 0;

		if(C_Sens2Str == "CALL" || C_Sens2Str == "C")
			C_type2sens = 1;
		else if(C_Sens2Str == "PUT" || C_Sens2Str == "P")
			C_type2sens = -1;
		else
			C_type2sens = 0;

		double C_NbSimul,C_AsOf;
		XL_readNumCell(XL_AsOf, C_AsOf, "ARM_ERR: asof numeric expected",C_result);
		XL_readNumCell(XL_NbSimul, C_NbSimul, "ARM_ERR: asof numeric expected",C_result);
		int C_intNbSimul = (int)C_NbSimul;
		
		CCString C_SabrStr;
		XL_readStrCellWD(XL_SabrDiff,C_SabrStr,"Y","ARM_ERR : sabr diff : string expected", C_result);
		C_SabrStr.toUpper();
		bool sabrdiff = C_SabrStr == "Y" ? true : false;

		double C_CorrVol,C_CorrMeanRev, def = 0.;
		XL_readNumCellWD(XL_CorrMeanRev ,C_CorrMeanRev,def,"ARM_ERR : sabr diff : string expected", C_result);
		XL_readNumCellWD(XL_CorrVol,C_CorrVol,def,"ARM_ERR : sabr diff : string expected", C_result);

		exportFunc26Args<	double,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							int,
							bool,
							vector<double>,
							vector<double>,
							int,
							int,
							int,
							double,
							double,
							vector<double>,
							vector<double>
						> 
			ourFunc(
				C_AsOf,
				C_ResetDates,
				C_Fwds1,
				C_FwdsVols1,
				C_FwdsPartVols1,
				C_Vols1,
				C_VVols1,
				C_Rho1,
				C_Fwds2,
				C_FwdsVols2,
				C_FwdsPartVols2,
				C_Vols2,
				C_VVols2,
				C_Rho2,
				C_R12Correl,
				C_intNbSimul,
				sabrdiff,
				C_fwdLev,
				C_fwdStrikes,
				C_typevb,
				C_type1sens,
				C_type2sens,
				C_CorrMeanRev,
				C_CorrVol,
				C_Levier,
				C_Fixed,
				ARMLOCAL_SpreadVBProxy_Create);

		/// call the general function
		fillXL_Result( LOCAL_PROXY_CLASS, ourFunc, C_result, XL_result, true );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_SpreadVBProxy_GetInfo(
	LPXLOPER XL_SpreadVBProxyId,
	LPXLOPER XL_Info
	)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		CCString C_ProxyStrId;
		XL_readStrCell( XL_SpreadVBProxyId, C_ProxyStrId,	" ARM_ERR: Proxy Id: Object expected", C_result);
		long C_ProxyId = LocalGetNumObjectId(C_ProxyStrId);
		
		CCString C_InfoStr;
		XL_readStrCell( XL_Info, C_InfoStr, "ARM_ERR : Info : string expected", C_result);
		
		int info;
		C_InfoStr.toUpper();

		if(C_InfoStr == "VB1")
			info = 1;
		else if(C_InfoStr == "VB2")
			info = 2;
		else if(C_InfoStr == "COUPON")
			info = 3;
		else
			info = 0;

		vector<double> vresult;

		long retCode = ARMLOCAL_SpreadVBProxy_GetInfo(C_ProxyId,info,vresult,C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
			XL_writeNumVector(XL_result, vresult, " ARM_ERR: Could not get result data", C_result);
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_Berm2DatesProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_strike,
	LPXLOPER XL_DFs,
	LPXLOPER XL_fwdRatesVol,
	LPXLOPER XL_fwdRatesPartVol,
	LPXLOPER XL_vols,
	LPXLOPER XL_volvols,
	LPXLOPER XL_rho,
	LPXLOPER XL_Rho12,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_TypeDiff)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		vector<double> defaultVector(0);
		defaultVector.clear();

		vector<double> C_ResetDates;
		XL_readNumVector(XL_resetDates,C_ResetDates," ARM_ERR: ResetDates: array of numeric expected",C_result);
		vector<double> C_Fwds;
		XL_readNumVector(XL_fwdRates,C_Fwds," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_DFs;
		XL_readNumVector(XL_DFs,C_DFs," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_FwdsVols;
		XL_readNumVector(XL_fwdRatesVol,C_FwdsVols," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_FwdsPartVols;
		XL_readNumVector(XL_fwdRatesPartVol,C_FwdsPartVols," ARM_ERR: FwdRates: array of numeric expected",C_result);
		vector<double> C_Vols;
		XL_readNumVector(XL_vols,C_Vols," ARM_ERR: Vols: array of numeric expected",C_result);
		vector<double> C_VVols;
		XL_readNumVector(XL_volvols,C_VVols," ARM_ERR: Vols: array of numeric expected",C_result);
		vector<double> C_Rho;
		XL_readNumVector(XL_rho,C_Rho," ARM_ERR: Vols: array of numeric expected",C_result);
		double C_Strike;
		XL_readNumCell(XL_strike, C_Strike, " ARM_ERR : Strike numeric expected", C_result);
		double C_Rho12;
		XL_readNumCell(XL_Rho12, C_Rho12, " ARM_ERR : Rate Correl numeric expected", C_result);

		double C_NbSimul,C_AsOf;
		XL_readNumCell(XL_AsOf, C_AsOf, "ARM_ERR: asof numeric expected",C_result);
		XL_readNumCell(XL_nbSimul, C_NbSimul, "ARM_ERR: asof numeric expected",C_result);
		int C_intNbSimul = (int)C_NbSimul;
		
		CCString C_DiffStr;
		XL_readStrCell(XL_TypeDiff,C_DiffStr,"ARM_ERR : sabr diff : string expected", C_result);
		C_DiffStr.toUpper();
		int C_DiffType;
		if(C_DiffStr == "SABR2")
			C_DiffType = 3;
		else if(C_DiffStr == "SABR")
			C_DiffType = 2;
		else if(C_DiffStr == "BLACK")
			C_DiffType = 1;
		else
			C_DiffType = 0;

		exportFunc13Args<	double,
							vector<double>,
							vector<double>,
							double,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							vector<double>,
							double,
							int,
							int
						> 
			ourFunc(
				C_AsOf,
				C_ResetDates,
				C_Fwds,
				C_Strike,
				C_DFs,
				C_FwdsVols,
				C_FwdsPartVols,
				C_Vols,
				C_VVols,
				C_Rho,
				C_Rho12,
				C_intNbSimul,
				C_DiffType,
				ARMLOCAL_Berm2DatesProxy_Create);

		/// call the general function
		fillXL_Result( LOCAL_PROXY_CLASS, ourFunc, C_result, XL_result, true );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Local_BiSABR_Digital_SpreadOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Digital_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_Digital_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	


	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_Digital_SpreadOption(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_K,
		C_T,
		C_CallPut,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_Digital_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Local_BiSABR_DigitalPaysS1_SpreadOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_DigitalPaysS1_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_DigitalPaysS1_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	


	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_DigitalPaysS1_SpreadOption(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_K,
		C_T,
		C_CallPut,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_DigitalPaysS1_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Local_BiSABR_DigitalPaysS2_SpreadOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_DigitalPaysS2_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_DigitalPaysS2_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	


	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_DigitalPaysS2_SpreadOption(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_K,
		C_T,
		C_CallPut,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_DigitalPaysS2_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Local_BiSABR_DigitalPaysS3_SpreadOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_DigitalPaysS3_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_S3params,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_DigitalPaysS3_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	
	vector<double> C_S3params;
	XL_readNumVector(XL_S3params,C_S3params," ARM_ERR: initialparams Vector: array of numeric expected",C_result);


	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_DigitalPaysS3_SpreadOption(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_S3params,
		C_K,
		C_T,
		C_CallPut,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_DigitalPaysS3_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Quantile and distribution
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Quantile(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	


	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_Quantile(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_K,
		C_T,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_Quantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Distribution(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_flag
	)
		{
			ADD_LOG("Local_BiSABR_Distribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,				" ARM_ERR: F1 numeric expected",C_result);
	double C_alpha1;
	XL_readNumCell(XL_alpha1,C_alpha1,				" ARM_ERR: alpha1 numeric expected",C_result);	
	double C_beta1;
	XL_readNumCell(XL_beta1,C_beta1,				" ARM_ERR: beta1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,				" ARM_ERR: nu1 numeric expected",C_result);	
	double C_F2;
	XL_readNumCell(XL_F2,C_F2,				" ARM_ERR: F2 numeric expected",C_result);	
	double C_alpha2;
	XL_readNumCell(XL_alpha2,C_alpha2,			" ARM_ERR: alpha2 numeric expected",C_result);	
	double C_beta2;
	XL_readNumCell(XL_beta2,C_beta2,			" ARM_ERR: beta2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhov numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_BiSABR_Distribution(
		C_F1,
		C_alpha1,
		C_beta1,
		C_rho1,
		C_nu1,
		C_F2,
		C_alpha2,
		C_beta2,
		C_rho2,
		C_nu2,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_K,
		C_T,
		C_flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_Distribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BetaEqualZeroSABR(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_mu,
	LPXLOPER XL_alpha,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_CallPut
	)
		{
			ADD_LOG("Local_BetaEqualZeroSABR");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F;
	XL_readNumCell(XL_F,C_F,				" ARM_ERR: F numeric expected",C_result);
	double C_K;
	XL_readNumCell(XL_K,C_K,				" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,				" ARM_ERR: T numeric expected",C_result);	
	double C_mu;
	XL_readNumCell(XL_mu,C_mu,				" ARM_ERR: mu numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,		" ARM_ERR: alpha numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,			" ARM_ERR: rho numeric expected",C_result);	
	double C_nu;
	XL_readNumCell(XL_nu,C_nu,				" ARM_ERR: nu numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);



	/// call the function
	long retCode=ARMLOCAL_BetaEqualZeroSABR(
		C_F,
		C_K,
		C_T,
		C_mu,
		C_alpha,
		C_rho,
		C_nu,
		C_CallPut,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BetaEqualZeroSABR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Shifted2LogNormal_Distribution(
	LPXLOPER XL_f1,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_f2,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha,
	LPXLOPER XL_rho,
	LPXLOPER XL_T,
	LPXLOPER XL_n,
	LPXLOPER XL_K
	)
		{
			ADD_LOG("Local_Shifted2LogNormal_Distribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,				" ARM_ERR: f1 numeric expected",C_result);
	double C_sigma1;
	XL_readNumCell(XL_sigma1,C_sigma1,				" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,				" ARM_ERR: f2 numeric expected",C_result);	
	double C_sigma2;
	XL_readNumCell(XL_sigma2,C_sigma2,				" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,				" ARM_ERR: alpha numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,				" ARM_ERR: rho numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	
	double C_n;
	XL_readNumCell(XL_n,C_n,			" ARM_ERR: n numeric expected",C_result);	

	/// call the function
	long retCode=ARMLOCAL_Shifted2LogNormal_Distribution(
		C_f1,
		C_sigma1,
		C_f2,
		C_sigma2,
		C_alpha,
		C_rho,
		C_K,
		C_T,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Shifted2LogNormal_Distribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Shifted2LogNormal_Quantile(
	LPXLOPER XL_f1,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_f2,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha,
	LPXLOPER XL_rho,
	LPXLOPER XL_T,
	LPXLOPER XL_n,
	LPXLOPER XL_K
	)
		{
			ADD_LOG("Local_Shifted2LogNormal_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_f1;
	XL_readNumCell(XL_f1,C_f1,				" ARM_ERR: f1 numeric expected",C_result);
	double C_sigma1;
	XL_readNumCell(XL_sigma1,C_sigma1,				" ARM_ERR: sigma1 numeric expected",C_result);	
	double C_f2;
	XL_readNumCell(XL_f2,C_f2,				" ARM_ERR: f2 numeric expected",C_result);	
	double C_sigma2;
	XL_readNumCell(XL_sigma2,C_sigma2,				" ARM_ERR: sigma2 numeric expected",C_result);	
	double C_alpha;
	XL_readNumCell(XL_alpha,C_alpha,				" ARM_ERR: alpha numeric expected",C_result);	
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,				" ARM_ERR: rho numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,			" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	
	double C_n;
	XL_readNumCell(XL_n,C_n,			" ARM_ERR: n numeric expected",C_result);	

	/// call the function
	long retCode=ARMLOCAL_Shifted2LogNormal_Quantile(
		C_f1,
		C_sigma1,
		C_f2,
		C_sigma2,
		C_alpha,
		C_rho,
		C_K,
		C_T,
		C_n,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Shifted2LogNormal_Quantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Local_BiSABR_S3_SpreadOption  pays A1*(S1-S2) +B1*S3 -K1 if A2*(S1-S2) +B2*S3 -K2 >0
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_S3_SpreadOption(
	LPXLOPER XL_S1params,
	LPXLOPER XL_S2params,
	LPXLOPER XL_S3params,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_Correlation,
	LPXLOPER XL_T,
	LPXLOPER XL_A1,
	LPXLOPER XL_B1,
	LPXLOPER XL_K1,
	LPXLOPER XL_A2,
	LPXLOPER XL_B2,
	LPXLOPER XL_K2,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps
	)
		{
			ADD_LOG("Local_BiSABR_S3_SpreadOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_S1params;
	XL_readNumVector(XL_S1params,C_S1params," ARM_ERR: initialparams Vector S1: array of numeric expected",C_result);
	vector<double> C_S2params;
	XL_readNumVector(XL_S2params,C_S2params," ARM_ERR: initialparams Vector S2: array of numeric expected",C_result);
	vector<double> C_S3params;
	XL_readNumVector(XL_S3params,C_S3params," ARM_ERR: initialparams Vector S3: array of numeric expected",C_result);


	double C_rhos;
	XL_readNumCell(XL_rhos,C_rhos,			" ARM_ERR: rhos numeric expected",C_result);	
	double C_rhov;
	XL_readNumCell(XL_rhov,C_rhov,			" ARM_ERR: rhovnumeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_Correlation;
	XL_readNumCell(XL_Correlation,C_Correlation,			" ARM_ERR: Correlation numeric expected",C_result);	


	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	

	double C_A1;
	XL_readNumCell(XL_A1,C_A1,		" ARM_ERR: A1 numeric expected",C_result);	
	double C_B1;
	XL_readNumCell(XL_B1,C_B1,		" ARM_ERR: B1 numeric expected",C_result);	
	double C_K1;
	XL_readNumCell(XL_K1,C_K1,		" ARM_ERR: K1numeric expected",C_result);	
	double C_A2;
	XL_readNumCell(XL_A2,C_A2,		" ARM_ERR: A2 numeric expected",C_result);	
	double C_B2;
	XL_readNumCell(XL_B2,C_B2,		" ARM_ERR: B2 numeric expected",C_result);	
	double C_K2;
	XL_readNumCell(XL_K2,C_K2,		" ARM_ERR: K2 numeric expected",C_result);	

	double C_flag;
	XL_readNumCell(XL_flag,C_flag,			" ARM_ERR: flag numeric expected",C_result);
	double C_nbsteps;
	XL_readNumCell(XL_nbsteps,C_nbsteps,			" ARM_ERR: nbsteps numeric expected",C_result);	



	/// call the function
	long retCode=ARMLOCAL_BiSABR_S3_SpreadOption(
		C_S1params,
		C_S2params,
		C_S3params,
		C_rhos,
		C_rhov,
		C_rhoc12,
		C_rhoc21,
		C_Correlation,
		C_T,
		C_A1,
		C_B1,
		C_K1,
		C_A2,
		C_B2,
		C_K2,
		C_flag,
		C_nbsteps,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSABR_S3_SpreadOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Heston_OptionPrice(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_InitialVar,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Times,
	LPXLOPER XL_Levels
	)
{
	ADD_LOG("Local_Heston_OptionPrice");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_AsOfDate;
		XL_readNumCell(XL_AsOfDate, C_AsOfDate, "ARM_ERR : asofdate : numeric expeted", C_result);
		double C_ResetTime;
		XL_readNumCell(XL_ResetTime, C_ResetTime, "ARM_ERR : reset time: numeric expeted", C_result);
		double C_Forward;
		XL_readNumCell(XL_Forward, C_Forward, "ARM_ERR : forward: numeric expeted", C_result);
		double C_Strike;
		XL_readNumCell(XL_Strike,  C_Strike, "ARM_ERR : strike: numeric expeted", C_result);
		double C_CallPut;
		XL_GETCONVCALLORPUT(XL_CallOrPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
		double C_InitialVar;
		XL_readNumCell(XL_InitialVar, C_InitialVar, "ARM_ERR : initial var: numeric expeted", C_result);
		double C_Kappa;
		XL_readNumCell(XL_Kappa, C_Kappa , "ARM_ERR : kappa : numeric expeted", C_result);
		double C_Rho;
		XL_readNumCell(XL_Rho, C_Rho, "ARM_ERR : rho : numeric expeted", C_result);
		double C_Theta;
		XL_readNumCell(XL_Theta, C_Theta, "ARM_ERR : theta : numeric or array expeted", C_result);
		double C_VVol;
		XL_readNumCell(XL_VVol, C_VVol, "ARM_ERR : vvol : numeric or array expeted", C_result);
		vector<double> C_DefaultVect, C_Times;
		XL_readNumVectorWD(XL_Times, C_Times, C_DefaultVect, "ARM_ERR : times : numeric or array expeted", C_result);
		double C_Shift, C_DefaultShift = 1.;
		XL_readNumCellWD(XL_Shift, C_Shift, C_DefaultShift, "ARM_ERR : shift : numeric expeted", C_result);
		vector<double> C_Level,C_DefaultLevel(1,1.);
		XL_readNumVectorWD(XL_Levels, C_Level, C_DefaultLevel, "ARM_ERR : level : numeric or array expeted", C_result);
		
		long retCode=ARMLOCAL_Heston_OptionPrice(C_AsOfDate, C_ResetTime, C_Forward, C_Strike, (int)C_CallPut, C_InitialVar, C_Kappa, 
			C_Rho, C_Theta, C_VVol, C_Shift, C_Times, C_Level, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Heston2B_OptionPrice(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_InitialVar1,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Theta1,
	LPXLOPER XL_VVol1,
	LPXLOPER XL_InitialVar2,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Theta2,
	LPXLOPER XL_VVol2,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Times,
	LPXLOPER XL_Levels
	)
{
	ADD_LOG("Local_Heston2B_OptionPrice");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_AsOfDate;
		XL_readNumCell(XL_AsOfDate, C_AsOfDate, "ARM_ERR : asofdate : numeric expeted", C_result);
		double C_ResetTime;
		XL_readNumCell(XL_ResetTime, C_ResetTime, "ARM_ERR : reset time: numeric expeted", C_result);
		double C_Forward;
		XL_readNumCell(XL_Forward, C_Forward, "ARM_ERR : forward: numeric expeted", C_result);
		double C_Strike;
		XL_readNumCell(XL_Strike,  C_Strike, "ARM_ERR : strike: numeric expeted", C_result);
		double C_CallPut;
		XL_GETCONVCALLORPUT(XL_CallOrPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
		double C_InitialVar1;
		XL_readNumCell(XL_InitialVar1, C_InitialVar1, "ARM_ERR : initial var1: numeric expeted", C_result);
		double C_Kappa1;
		XL_readNumCell(XL_Kappa1, C_Kappa1 , "ARM_ERR : kappa1 : numeric expeted", C_result);
		double C_Rho1;
		XL_readNumCell(XL_Rho1, C_Rho1, "ARM_ERR : rho1 : numeric expeted", C_result);
		double C_Theta1;
		XL_readNumCell(XL_Theta1, C_Theta1, "ARM_ERR : theta1 : numeric or array expeted", C_result);
		double C_VVol1;
		XL_readNumCell(XL_VVol1, C_VVol1, "ARM_ERR : vvol1 : numeric or array expeted", C_result);
		double C_InitialVar2;
		XL_readNumCell(XL_InitialVar2, C_InitialVar2, "ARM_ERR : initial var2: numeric expeted", C_result);
		double C_Kappa2;
		XL_readNumCell(XL_Kappa2, C_Kappa2 , "ARM_ERR : kappa2 : numeric expeted", C_result);
		double C_Rho2;
		XL_readNumCell(XL_Rho2, C_Rho2, "ARM_ERR : rho2 : numeric expeted", C_result);
		double C_Theta2;
		XL_readNumCell(XL_Theta2, C_Theta2, "ARM_ERR : theta2 : numeric or array expeted", C_result);
		double C_VVol2;
		XL_readNumCell(XL_VVol2, C_VVol2, "ARM_ERR : vvol2 : numeric or array expeted", C_result);

		vector<double> C_DefaultVect, C_Times;
		XL_readNumVectorWD(XL_Times, C_Times, C_DefaultVect, "ARM_ERR : times : numeric or array expeted", C_result);
		double C_Shift, C_DefaultShift = 1.;
		XL_readNumCellWD(XL_Shift, C_Shift, C_DefaultShift, "ARM_ERR : shift : numeric expeted", C_result);
		vector<double> C_Level,C_DefaultLevel(1,1.);
		XL_readNumVectorWD(XL_Levels, C_Level, C_DefaultLevel, "ARM_ERR : level : numeric or array expeted", C_result);
		
		long retCode=ARMLOCAL_Heston2b_OptionPrice(C_AsOfDate, C_ResetTime, C_Forward, C_Strike, (int)C_CallPut, 
			C_InitialVar1, C_Kappa1, C_Rho1, C_Theta1, C_VVol1,
			C_InitialVar2, C_Kappa2, C_Rho2, C_Theta2, C_VVol2,
			C_Shift, C_Times, C_Level, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_MixteHeston_OptionPrice(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_Sigma,
	LPXLOPER XL_InitialVar,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Times,
	LPXLOPER XL_Levels
	)
{
	ADD_LOG("Local_MixteHeston_OptionPrice");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_AsOfDate;
		XL_readNumCell(XL_AsOfDate, C_AsOfDate, "ARM_ERR : asofdate : numeric expeted", C_result);
		double C_ResetTime;
		XL_readNumCell(XL_ResetTime, C_ResetTime, "ARM_ERR : reset time: numeric expeted", C_result);
		double C_Forward;
		XL_readNumCell(XL_Forward, C_Forward, "ARM_ERR : forward: numeric expeted", C_result);
		double C_Strike;
		XL_readNumCell(XL_Strike,  C_Strike, "ARM_ERR : strike: numeric expeted", C_result);
		double C_CallPut;
		XL_GETCONVCALLORPUT(XL_CallOrPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
		double C_Sigma;
		XL_readNumCell(XL_Sigma,  C_Sigma, "ARM_ERR : sigma: numeric expeted", C_result);
		double C_InitialVar;
		XL_readNumCell(XL_InitialVar, C_InitialVar, "ARM_ERR : initial var: numeric expeted", C_result);
		double C_Kappa;
		XL_readNumCell(XL_Kappa, C_Kappa , "ARM_ERR : kappa : numeric expeted", C_result);
		double C_Rho;
		XL_readNumCell(XL_Rho, C_Rho, "ARM_ERR : rho : numeric expeted", C_result);
		double C_Theta;
		XL_readNumCell(XL_Theta, C_Theta, "ARM_ERR : theta : numeric or array expeted", C_result);
		double C_VVol;
		XL_readNumCell(XL_VVol, C_VVol, "ARM_ERR : vvol : numeric or array expeted", C_result);
		vector<double> C_DefaultVect, C_Times;
		XL_readNumVectorWD(XL_Times, C_Times, C_DefaultVect, "ARM_ERR : times : numeric or array expeted", C_result);
		double C_Shift, C_DefaultShift = 1.;
		XL_readNumCellWD(XL_Shift, C_Shift, C_DefaultShift, "ARM_ERR : shift : numeric expeted", C_result);
		vector<double> C_Level,C_DefaultLevel(1,1.);
		XL_readNumVectorWD(XL_Levels, C_Level, C_DefaultLevel, "ARM_ERR : level : numeric or array expeted", C_result);
		
		long retCode=ARMLOCAL_MixteHeston_OptionPrice(C_AsOfDate, C_ResetTime, C_Forward, C_Strike, 
			(int)C_CallPut, C_Sigma, C_InitialVar, C_Kappa, 
			C_Rho, C_Theta, C_VVol, C_Shift, C_Times, C_Level, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	
///			Normal Heston Vanilla option
///
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Normal_Heston_VanillaCall(
	LPXLOPER XL_rho,
	LPXLOPER XL_lambdaV,
	LPXLOPER XL_thetaV,
	LPXLOPER XL_kappaV,
	LPXLOPER XL_V0,
	LPXLOPER XL_S0,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_lambdaB,
	LPXLOPER XL_callput,
	LPXLOPER XL_nbfirst,
	LPXLOPER XL_nb,
	LPXLOPER XL_NbStage,
	LPXLOPER XL_NbOscill,
	LPXLOPER XL_prec
	)
		{
			ADD_LOG("Local_Normal_Heston_VanillaCall");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_rho;
	XL_readNumCell(XL_rho,C_rho,				" ARM_ERR: rho numeric expected",C_result);
	double C_lambdaV;
	XL_readNumCell(XL_lambdaV,C_lambdaV,				" ARM_ERR: lambdaV numeric expected",C_result);	
	double C_thetaV;
	XL_readNumCell(XL_thetaV,C_thetaV,				" ARM_ERR: thetaV numeric expected",C_result);	
	double C_kappaV;
	XL_readNumCell(XL_kappaV,C_kappaV,				" ARM_ERR: kappaV numeric expected",C_result);	
	double C_V0;
	XL_readNumCell(XL_V0,C_V0,				" ARM_ERR: V0 numeric expected",C_result);	
	double C_S0;
	XL_readNumCell(XL_S0,C_S0,				" ARM_ERR: S0 numeric expected",C_result);	
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	
	double C_lambdaB;
	double C_lambdaB_def=0.1;
	XL_readNumCellWD(XL_lambdaB,C_lambdaB,C_lambdaB_def,			" ARM_ERR: lambdaB numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_nbfirst;
	double C_nbfirst_def=40;
	XL_readNumCellWD(XL_nbfirst,C_nbfirst,C_nbfirst_def,	" ARM_ERR: nbfirst numeric expected",C_result);	
	double C_nb;
	double C_nb_def=40;
	XL_readNumCellWD(XL_nb,C_nb,C_nb_def,			" ARM_ERR: nb numeric expected",C_result);
	double C_NbStage;
	double C_NbStage_def=-1; /// to force the automatique determination of stage
	XL_readNumCellWD(XL_NbStage,C_NbStage,C_NbStage_def,			" ARM_ERR: NbStage numeric expected",C_result);	
	double C_NbOscill;
	double C_NbOscill_def=0; // to force the default = 0.01
	XL_readNumCellWD(XL_NbOscill,C_NbOscill,C_NbOscill_def,			" ARM_ERR: NbOscill numeric expected",C_result);	
	double C_prec;
	double C_prec_def=1e-8;
	XL_readNumCellWD(XL_prec,C_prec,C_prec_def,			" ARM_ERR: prec numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Normal_Heston_VanillaCall(
		C_rho,
		C_lambdaV,
		C_thetaV,
		C_kappaV,
		C_V0,
		C_S0,
		C_k,
		C_T,
		C_lambdaB,
		C_callput,
		C_nbfirst,
		C_nb,
		C_NbStage,
		C_NbOscill,
		C_prec,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_Heston_VanillaCall" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_SABR_To_Heston_SmileCalibration_Create(
	LPXLOPER XL_resetTime,
	LPXLOPER XL_fwdRate,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_RhoSABR,
	LPXLOPER XL_Nu,
	LPXLOPER XL_Sabr_Type,
	LPXLOPER XL_InitialVar,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Shift)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_resetTime;
		XL_readNumCell(XL_resetTime, C_resetTime, "ARM_ERR : reset : numeric expeted", C_result);
		double C_fwdRate;
		XL_readNumCell(XL_fwdRate, C_fwdRate, "ARM_ERR : fwdrate: numeric expeted", C_result);
		double C_ATMVol;
		XL_readNumCell(XL_ATMVol, C_ATMVol, "ARM_ERR : ATMVol: numeric expeted", C_result);
		double C_Alpha;
		XL_readNumCell(XL_Alpha, C_Alpha, "ARM_ERR : Alpha: numeric expeted", C_result);
		double C_Beta;
		XL_readNumCell(XL_Beta, C_Beta, "ARM_ERR : Beta: numeric expeted", C_result);
		double C_RhoSABR;
		XL_readNumCell(XL_RhoSABR, C_RhoSABR, "ARM_ERR : RhoSABR: numeric expeted", C_result);
		double C_Nu;
		XL_readNumCell(XL_Nu, C_Nu, "ARM_ERR : Nu: numeric expeted", C_result);
		double C_Sabr_Type;
		XL_GETCONVSABRFLAG(XL_Sabr_Type,C_Sabr_Type," ARM_ERR: SABR flag string expected",C_result);
		
		double C_InitialVar;
		XL_readNumCell(XL_InitialVar, C_InitialVar, "ARM_ERR : InitialVar: numeric expeted", C_result);
		double C_Kappa;
		XL_readNumCell(XL_Kappa, C_Kappa, "ARM_ERR : Kappa: numeric expeted", C_result);
		double dbldef = -999.;
		double C_Theta;
		XL_readNumCellWD(XL_Theta, C_Theta, dbldef, "ARM_ERR : theta : numeric expeted", C_result);
		double C_Rho;
		XL_readNumCellWD(XL_Rho, C_Rho, dbldef, "ARM_ERR : Rho : numeric expeted", C_result);
		double C_Shift;
		XL_readNumCellWD(XL_Shift, C_Shift, dbldef, "ARM_ERR : Shift : numeric expeted", C_result);

		bool C_CalibTheta = C_Theta == dbldef;
		bool C_CalibRho = C_Rho == dbldef;
		bool C_CalibShift = C_Shift == dbldef;

		exportFunc16Args<	double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							bool,
							bool,
							bool
						> 
			ourFunc(
				C_resetTime,
				C_fwdRate,
				C_ATMVol,
				C_Alpha,
				C_Beta,
				C_RhoSABR,
				C_Nu,
				C_Sabr_Type,
				C_InitialVar,
				C_Kappa,
				C_Theta,
				C_Rho,
				C_Shift,
				C_CalibTheta,
				C_CalibRho,
				C_CalibShift,
				ARMLOCAL_SABR_To_Heston_SmileCalibration_Create);

		/// call the general function
		fillXL_Result( LOCAL_SABRHESTONCALIB_CLASS, ourFunc, C_result, XL_result, true );
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_SABR_To_Heston_SmileCalibration_GetValue(
	LPXLOPER XL_CalibrationId,
	LPXLOPER XL_Info
	)
{
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		CCString C_CalibStrId;
		XL_readStrCell( XL_CalibrationId, C_CalibStrId,	" ARM_ERR: Calib Id: Object expected", C_result);
		long C_CalibId = LocalGetNumObjectId(C_CalibStrId);
		
		CCString C_InfoStr;
		XL_readStrCell( XL_Info, C_InfoStr, "ARM_ERR : Info : string expected", C_result);
		
		int info;
		C_InfoStr.toUpper();

		if(C_InfoStr == "V0")
			info = 0;
		else if(C_InfoStr == "KAPPA")
			info = 1;
		else if(C_InfoStr == "THETA")
			info = 2;
		else if(C_InfoStr == "VVOL")
			info = 3;
		else if(C_InfoStr == "RHO")
			info = 4;
		else if(C_InfoStr == "SHIFT")
			info = 5;
		else if(C_InfoStr == "LEVEL")
			info = 6;
		else if(C_InfoStr == "ERROR")
			info = 7;

		long retCode = ARMLOCAL_SABR_To_Heston_SmileCalibration_GetValue(C_CalibId,info,C_result);

	    /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
	    }
	    else
	    {
		    ARM_ERR();
	    }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
_declspec(dllexport) LPXLOPER WINAPI Local_SmileCalibration(
	LPXLOPER	XL_AsOfDate,
	LPXLOPER	XL_CalibTimes,
	LPXLOPER	XL_Forwards,
	LPXLOPER	XL_MktVols,
	LPXLOPER	XL_Strikes,
	LPXLOPER	XL_CalibParamId,
	LPXLOPER	XL_ConstraintStrikes,
	LPXLOPER	XL_ConstraintVols,
	LPXLOPER	XL_Weights)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_AsOfDate;
		XL_readNumCell(XL_AsOfDate, C_AsOfDate, "ARM_ERR : asofdate : numeric expeted", C_result);
		vector<double> C_CalibTimes;
		XL_readNumVector(XL_CalibTimes, C_CalibTimes, "ARM_ERR : calib time : numeric or array expeted", C_result);
		vector<double> C_Forwards;
		XL_readNumVector(XL_Forwards, C_Forwards, "ARM_ERR : forward : numeric or array expeted", C_result);
		vector< vector<double> > C_MktVols;
		XL_readNumVectorVector(XL_MktVols, C_MktVols, "ARM_ERR : mkt vols : array of numeric expected", C_result);
		vector< vector<double> > C_Strikes;
		XL_readNumVectorVector(XL_Strikes, C_Strikes," ARM_ERR: strikes : array of numeric expected",C_result);

		vector<double> C_ConstraintVols, C_ConstraintStrikes, C_Weights, defVec(0);

		XL_readNumVectorWD(XL_ConstraintStrikes, C_ConstraintStrikes, defVec, " ARM_ERR: constraint strikes : numeric or array expected", C_result);
		XL_readNumVectorWD(XL_ConstraintVols, C_ConstraintVols, defVec, " ARM_ERR: constraint vols : numeric or array expected", C_result);
		XL_readNumVectorWD(XL_Weights, C_Weights, defVec, " ARM_ERR: constraint strikes : numeric or array expected", C_result);

		bool calibConstraint = C_ConstraintStrikes.size() == 0 || C_ConstraintVols.size() == 0 ? false : true;
		if(C_Weights.size() == 0) C_Weights.resize(C_CalibTimes.size(), 1.);

		CCString C_CalibStrId;
		XL_readStrCell( XL_CalibParamId, C_CalibStrId,	" ARM_ERR: Calib Param Id: Object expected", C_result);
		long C_CalibId = LocalGetNumObjectId(C_CalibStrId);

		vector<double> outResult;
		std::deque<bool> outResultBool;
		int rows, cols;

		long retCode=ARMLOCAL_SmileCalibration(
			C_AsOfDate,
			C_CalibTimes,
			C_Forwards,
			C_MktVols,
			C_Strikes,
			C_CalibId,
			C_ConstraintStrikes,
			C_ConstraintVols,
			calibConstraint,
			C_Weights,
			rows, cols, outResult, outResultBool,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			VECTOR<CCString> C_OutResult;
			C_OutResult.resize(rows * cols, "");
			VECTOR<long> types(C_OutResult.size(), xltypeStr);

			int i, k;
			char outstr[30];

			for(i = 0; i < rows; i++)
			{
				for(k = 0; k < cols; k++)
				{
					if(outResultBool[i+k*rows])
					{
						sprintf(outstr,"%lf",outResult[i+k*rows]);
						C_OutResult[k + i*cols] = outstr;
						types[k + i*cols] = xltypeNum;
					}
				}
			}

			XL_writeStrMatrixSizeAndType(XL_result, C_OutResult, types, rows, cols, " ARM_ERR: Could not set the num matrix", C_result);
		}
		else
		{
			ARM_ERR();
		}

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_SmileCalib_Spread2Heston(
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Fwd1,
	LPXLOPER XL_MktVols1,
	LPXLOPER XL_Strikes1,
	LPXLOPER XL_ConstrVol1,
	LPXLOPER XL_ConstrK1,
	LPXLOPER XL_Fwd2,
	LPXLOPER XL_MktVols2,
	LPXLOPER XL_Strikes2,
	LPXLOPER XL_ConstrVol2,
	LPXLOPER XL_ConstrK2,
	LPXLOPER XL_MktVolsSpread,
	LPXLOPER XL_StrikesSpread,
	LPXLOPER XL_ConstrVolSpread,
	LPXLOPER XL_ConstrKSpread,
	LPXLOPER XL_v0,
	LPXLOPER XL_kappa,
	LPXLOPER XL_theta
)
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		double C_ResetTime, C_Fwd1, C_Fwd2;
		double C_ConstrVol1, C_ConstrK1, C_ConstrVol2, C_ConstrK2;
		double C_ConstrVolSpread, C_ConstrKSpread;
		double C_V0, C_Kappa, C_Theta, def = -999.;
		vector <double> C_MktVols1, C_Strikes1, C_MktVols2, C_Strikes2, C_MktVolsSpread, C_StrikesSpread;

		XL_readNumCell(XL_ResetTime, C_ResetTime, "ARM_ERR : reset: numeric expeted", C_result);
		XL_readNumCell(XL_Fwd1, C_Fwd1, "ARM_ERR : fwd1: numeric expeted", C_result);
		XL_readNumCell(XL_Fwd2, C_Fwd2, "ARM_ERR : fwd2: numeric expeted", C_result);
		XL_readNumCell(XL_ConstrVol1, C_ConstrVol1, "ARM_ERR : constr vol1: numeric expeted", C_result);
		XL_readNumCell(XL_ConstrK1, C_ConstrK1, "ARM_ERR : constr K1: numeric expeted", C_result);
		XL_readNumCell(XL_ConstrVol2, C_ConstrVol2, "ARM_ERR : constr vol2: numeric expeted", C_result);
		XL_readNumCell(XL_ConstrK2, C_ConstrK2, "ARM_ERR : constr K2: numeric expeted", C_result);
		XL_readNumCell(XL_ConstrVolSpread, C_ConstrVolSpread, "ARM_ERR : constr vol spread: numeric expeted", C_result);
		XL_readNumCell(XL_ConstrKSpread, C_ConstrKSpread, "ARM_ERR : constr K spread: numeric expeted", C_result);
		XL_readNumCell(XL_v0, C_V0, "ARM_ERR : v0: numeric expeted", C_result);
		XL_readNumCell(XL_kappa, C_Kappa, "ARM_ERR : kappa: numeric expeted", C_result);
		XL_readNumCellWD(XL_theta, C_Theta, def, "ARM_ERR : kappa: numeric expeted", C_result);
		bool C_calibTheta = C_Theta < K_DOUBLE_TOL ? true : false;

		XL_readNumVector(XL_MktVols1 , C_MktVols1 , "ARM_ERR : mkt vols1: numeric or array expeted", C_result);
		XL_readNumVector(XL_Strikes1 , C_Strikes1 , "ARM_ERR : strikes1: numeric or array expeted", C_result);
		XL_readNumVector(XL_MktVols2 , C_MktVols2 , "ARM_ERR : mkt vols2: numeric or array expeted", C_result);
		XL_readNumVector(XL_Strikes2 , C_Strikes2 , "ARM_ERR : strikes2: numeric or array expeted", C_result);
		XL_readNumVector(XL_MktVolsSpread , C_MktVolsSpread , "ARM_ERR : mkt vols spread: numeric or array expeted", C_result);
		XL_readNumVector(XL_StrikesSpread , C_StrikesSpread , "ARM_ERR : strikes spread: numeric or array expeted", C_result);

		vector<double> outResult;
		std::deque<bool> outResultBool;
		int rows, cols;

		long retCode=ARMLOCAL_SpreadSmileCalibration2Heston(
			C_ResetTime,
			C_Fwd1,
			C_MktVols1,
			C_Strikes1,
			C_ConstrVol1,
			C_ConstrK1,
			C_Fwd2,
			C_MktVols2,
			C_Strikes2,
			C_ConstrVol2,
			C_ConstrK2,
			C_MktVolsSpread,
			C_StrikesSpread,
			C_ConstrVolSpread,
			C_ConstrKSpread,
			C_V0,
			C_Kappa,
			C_Theta,
			C_calibTheta,
			rows, cols, outResult, outResultBool,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			VECTOR<CCString> C_OutResult;
			C_OutResult.resize(rows * cols, "");
			VECTOR<long> types(C_OutResult.size(), xltypeStr);

			int i, k;
			char outstr[30];

			for(i = 0; i < rows; i++)
			{
				for(k = 0; k < cols; k++)
				{
					if(outResultBool[i+k*rows])
					{
						sprintf(outstr,"%lf",outResult[i+k*rows]);
						C_OutResult[k + i*cols] = outstr;
						types[k + i*cols] = xltypeNum;
					}
				}
			}

			XL_writeStrMatrixSizeAndType(XL_result, C_OutResult, types, rows, cols, " ARM_ERR: Could not set the num matrix", C_result);
		}
		else
		{
			ARM_ERR();
		}

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_Spread2HestonVanilla(
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Fwd1,
	LPXLOPER XL_Fwd2,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_Nu,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Shift1,
	LPXLOPER XL_Shift2,
	LPXLOPER XL_Level1,
	LPXLOPER XL_Level2,
	LPXLOPER XL_Correl,
	LPXLOPER XL_Index1Lev,
	LPXLOPER XL_Index2Lev
)
{
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN

	{/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		double reset, fwd1, fwd2, strike, v0, theta, kappa, nu, rho1, rho2, shift1, shift2, level1, level2, correl;
		double index1lev, index2lev, def = 1.;

		XL_readNumCell(XL_ResetTime, reset, "ARM_ERR : reset: numeric expeted", C_result);
		XL_readNumCell(XL_Fwd1, fwd1, "ARM_ERR : fwd1: numeric expeted", C_result);
		XL_readNumCell(XL_Fwd2, fwd2, "ARM_ERR : fwd2: numeric expeted", C_result);
		XL_readNumCell(XL_Strike, strike, "ARM_ERR : fwd2: numeric expeted", C_result);
		XL_readNumCell(XL_V0, v0, "ARM_ERR : v0: numeric expeted", C_result);
		XL_readNumCell(XL_Kappa, kappa, "ARM_ERR : kappa: numeric expeted", C_result);
		XL_readNumCell(XL_Theta, theta, "ARM_ERR : theta : numeric expeted", C_result);
		XL_readNumCell(XL_Nu, nu, "ARM_ERR : theta : numeric expeted", C_result);
		XL_readNumCell(XL_Rho1, rho1, "ARM_ERR : rho1 : numeric expeted", C_result);
		XL_readNumCell(XL_Shift1, shift1, "ARM_ERR : shift1: numeric expeted", C_result);
		XL_readNumCell(XL_Level1, level1, "ARM_ERR : level1 : numeric expeted", C_result);
		XL_readNumCell(XL_Rho2, rho2, "ARM_ERR : rho2 : numeric expeted", C_result);
		XL_readNumCell(XL_Shift2, shift2, "ARM_ERR : shift2: numeric expeted", C_result);
		XL_readNumCell(XL_Level2, level2, "ARM_ERR : level2 : numeric expeted", C_result);
		XL_readNumCell(XL_Correl, correl, "ARM_ERR : correl : numeric expeted", C_result);
		XL_readNumCellWD(XL_Index1Lev, index1lev, def, "ARM_ERR : index1 leverage: numeric expeted", C_result);
		XL_readNumCellWD(XL_Index2Lev, index2lev, def, "ARM_ERR : index2 leverage: numeric expeted", C_result);

		double callPut;
		XL_GETCONVCALLORPUT(XL_CallPut,callPut," ARM_ERR: cap or floor string expected",C_result);

		long retCode=ARMLOCAL_Spread2HestonVanilla(reset, fwd1, fwd2, strike, callPut, v0,
						kappa, theta, nu, rho1, rho2, shift1, shift2, level1, level2, correl,
						index1lev, index2lev,
						C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TRiSABR_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_Spread2Heston_TOTEMCalibration(
	LPXLOPER XL_TOTEMMat,
	LPXLOPER XL_TOTEMStrikes,
	LPXLOPER XL_TOTEMPrices,
	LPXLOPER XL_FullScheduleReset,
	LPXLOPER XL_FullScheduleAnnuity,
	LPXLOPER XL_FullScheduleFwd1,
	LPXLOPER XL_FullScheduleFwd2,
	LPXLOPER XL_FwdCalibReset,
	LPXLOPER XL_LongFwds,
	LPXLOPER XL_LongVolsId,
	LPXLOPER XL_LongStrikesId,
	LPXLOPER XL_ShortFwds,
	LPXLOPER XL_ShortVolsId,
	LPXLOPER XL_ShortStrikesId,
	LPXLOPER XL_v0,
	LPXLOPER XL_kappa,
	LPXLOPER XL_theta,
	LPXLOPER XL_ConstrCorrel
)
{
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN

	{/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		vector<double> totemmat, totemstrikes, totemprices;
		vector<double> schedreset, schedannuity, schedfwd1, schedfwd2;
		vector<double> fwdcalibreset, longfwds, shortfwds;

		XL_readNumVector(XL_TOTEMMat, totemmat, "ARM_ERR : Totem Mat array expected", C_result);
		XL_readNumVector(XL_TOTEMStrikes, totemstrikes, "ARM_ERR : Totem Strikes array expected", C_result);
		XL_readNumVector(XL_TOTEMPrices, totemprices, "ARM_ERR : Totem Prices array expected", C_result);
		XL_readNumVector(XL_FullScheduleReset, schedreset, "ARM_ERR : Full Schedule Reset : array expected", C_result);
		XL_readNumVector(XL_FullScheduleAnnuity, schedannuity, "ARM_ERR : Full Schedule Annuity : array expected", C_result);
		XL_readNumVector(XL_FullScheduleFwd1, schedfwd1, "ARM_ERR : Full Schedule Spread : array expected", C_result);
		XL_readNumVector(XL_FullScheduleFwd2, schedfwd2, "ARM_ERR : Full Schedule Spread : array expected", C_result);
		XL_readNumVector(XL_FwdCalibReset, fwdcalibreset, "ARM_ERR : Fwds Calib Reset : array expected", C_result);
		XL_readNumVector(XL_LongFwds, longfwds, "ARM_ERR : Long Fwds : array expected", C_result);
		XL_readNumVector(XL_ShortFwds, shortfwds, "ARM_ERR : Short Fwds : array expected", C_result);

		double v0, kappa, theta, def = 0.;
		XL_readNumCell(XL_v0, v0, "ARM_ERR : v0 numeric expected", C_result);
		XL_readNumCell(XL_kappa, kappa, "ARM_ERR : kappa numeric expected", C_result);
		XL_readNumCellWD(XL_theta, theta, def, "ARM_ERR : theta numeric expected", C_result);

		CCString strConstr, defstr = "N";
		XL_readStrCellWD(XL_ConstrCorrel, strConstr, defstr, "ARM_ERR : Constrained Correl string (Y/N) expected", C_result);
		strConstr.toUpper();
		bool constrCorrel = strConstr == "N" ? false : true;

		long longVolId, longKId, shortVolId, shortKId;
		CCString strObjId;
		XL_readStrCell( XL_LongVolsId, strObjId, " ARM_ERR: Long Vols : Object (matrix) expected", C_result);
		longVolId = LocalGetNumObjectId(strObjId);
		XL_readStrCell( XL_LongStrikesId, strObjId, " ARM_ERR: Long Strikes : Object (matrix) expected", C_result);
		longKId = LocalGetNumObjectId(strObjId);
		XL_readStrCell( XL_ShortVolsId, strObjId, " ARM_ERR: Short Vols : Object (matrix) expected", C_result);
		shortVolId = LocalGetNumObjectId(strObjId);
		XL_readStrCell( XL_ShortStrikesId, strObjId, " ARM_ERR: Short Strikes : Object (matrix) expected", C_result);
		shortKId = LocalGetNumObjectId(strObjId);

		vector<double> outResult;
		std::deque<bool> outResultBool;
		int outRows, outCols;

		long retCode = ARMLOCAL_Spread2Heston_TOTEMCalibration(
			totemmat, 
			totemstrikes,
			totemprices,
			schedreset,
			schedannuity,
			schedfwd1,
			schedfwd2,
			fwdcalibreset,
			longfwds,
			longVolId,
			longKId,
			shortfwds,
			shortVolId,
			shortKId,
			v0,
			kappa,
			theta,
			constrCorrel,
			outRows, outCols, outResult, outResultBool, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			VECTOR<CCString> C_OutResult;
			C_OutResult.resize(outRows * outCols, "");
			VECTOR<long> types(C_OutResult.size(), xltypeStr);

			int i, k;
			char outstr[30];

			for(i = 0; i < outRows; i++)
			{
				for(k = 0; k < outCols; k++)
				{
					if(outResultBool[i+k*outRows] == true)
					{
						sprintf(outstr,"%lf",outResult[i+k*outRows]);
						C_OutResult[k + i*outCols] = outstr;
						types[k + i*outCols] = xltypeNum;
					}
					else
					{
						C_OutResult[k + i*outCols] = "";
						types[k + i*outCols] = xltypeStr;
					}
				}
			}

			XL_writeStrMatrixSizeAndType(XL_result, C_OutResult, types, outRows, outCols, " ARM_ERR: Could not set the num matrix", C_result);
		}
		else
		{
			ARM_ERR();
		}

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Calibrate" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	TriSABR : vanilla option 
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_TRiSABR_VanillaOption(
	LPXLOPER XL_S1params,
	LPXLOPER XL_S2params,
	LPXLOPER XL_S3params,
	LPXLOPER XL_rhos12,
	LPXLOPER XL_rhos23,
	LPXLOPER XL_rhos13,
	LPXLOPER XL_rhov12,
	LPXLOPER XL_rhov23,
	LPXLOPER XL_rhov13,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_rhoc23,
	LPXLOPER XL_rhoc32,
	LPXLOPER XL_rhoc13,
	LPXLOPER XL_rhoc31,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_callput,
	LPXLOPER XL_Sabr_Type
	)
		{
			ADD_LOG("Local_TRiSABR_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_S1params;
	XL_readNumVector(XL_S1params,C_S1params," ARM_ERR: initialparams Vector S1: array of numeric expected",C_result);
	vector<double> C_S2params;
	XL_readNumVector(XL_S2params,C_S2params," ARM_ERR: initialparams Vector S2: array of numeric expected",C_result);
	vector<double> C_S3params;
	XL_readNumVector(XL_S3params,C_S3params," ARM_ERR: initialparams Vector S3: array of numeric expected",C_result);
	double C_rhos12;
	XL_readNumCell(XL_rhos12,C_rhos12,			" ARM_ERR: rhos12 numeric expected",C_result);	
		double C_rhos23;
	XL_readNumCell(XL_rhos23,C_rhos23,			" ARM_ERR: rhos23 numeric expected",C_result);	
		double C_rhos13;
	XL_readNumCell(XL_rhos13,C_rhos13,			" ARM_ERR: rhos13 numeric expected",C_result);	
	double C_rhov12;
	XL_readNumCell(XL_rhov12,C_rhov12,			" ARM_ERR: rhov12 numeric expected",C_result);	
	double C_rhov23;
	XL_readNumCell(XL_rhov23,C_rhov23,			" ARM_ERR: rhov23 numeric expected",C_result);	
	double C_rhov13;
	XL_readNumCell(XL_rhov13,C_rhov13,			" ARM_ERR: rhov13 numeric expected",C_result);	
	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,			" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,			" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_rhoc23;
	XL_readNumCell(XL_rhoc23,C_rhoc23,			" ARM_ERR: rhoc23 numeric expected",C_result);	
	double C_rhoc32;
	XL_readNumCell(XL_rhoc32,C_rhoc32,			" ARM_ERR: rhoc32 numeric expected",C_result);	
	double C_rhoc13;
	XL_readNumCell(XL_rhoc13,C_rhoc13,			" ARM_ERR: rhoc13 numeric expected",C_result);	
	double C_rhoc31;
	XL_readNumCell(XL_rhoc31,C_rhoc31,			" ARM_ERR: rhoc31 numeric expected",C_result);	
	double C_K;
	XL_readNumCell(XL_K,C_K,		" ARM_ERR: K numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,		" ARM_ERR: T numeric expected",C_result);	
	double C_callput;
	XL_GETCONVCALLORPUT(XL_callput,C_callput," ARM_ERR: call or put string expected",C_result);
	double C_Sabr_Type;
	XL_GETCONVSABRFLAG(XL_Sabr_Type,C_Sabr_Type," ARM_ERR: SABR flag string expected",C_result);


	/// call the function
	long retCode=ARMLOCAL_TRiSABR_VanillaOption(
		C_S1params,
		C_S2params,
		C_S3params,
		C_rhos12,
		C_rhos23,
		C_rhos13,
		C_rhov12,
		C_rhov23,
		C_rhov13,
		C_rhoc12,
		C_rhoc21,
		C_rhoc23,
		C_rhoc32,
		C_rhoc13,
		C_rhoc31,
		C_K,
		C_T,
		C_callput,
		C_Sabr_Type,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TRiSABR_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_Util_TriSABR_Eigenvalues(
			LPXLOPER XL_rho1,   LPXLOPER XL_rho2,   LPXLOPER XL_rho3,
			LPXLOPER XL_rhos12, LPXLOPER XL_rhos23, LPXLOPER XL_rhos13,
			LPXLOPER XL_rhov12, LPXLOPER XL_rhov23, LPXLOPER XL_rhov13,
			LPXLOPER XL_rhoc12, LPXLOPER XL_rhoc13,
			LPXLOPER XL_rhoc21, LPXLOPER XL_rhoc23,
			LPXLOPER XL_rhoc31, LPXLOPER XL_rhoc32
	)
{
	ADD_LOG("Local_Util_TriSABR_Eigenvalues");
	

	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,				" ARM_ERR: rho1 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,				" ARM_ERR: rho2 numeric expected",C_result);	
	double C_rho3;
	XL_readNumCell(XL_rho3,C_rho3,				" ARM_ERR: rho3 numeric expected",C_result);	
	double C_rhos12;
	XL_readNumCell(XL_rhos12,C_rhos12,				" ARM_ERR: rhos12 numeric expected",C_result);	
	double C_rhos13;
	XL_readNumCell(XL_rhos13,C_rhos13,				" ARM_ERR: rhos13 numeric expected",C_result);	
	double C_rhos23;
	XL_readNumCell(XL_rhos23,C_rhos23,				" ARM_ERR: rhos23 numeric expected",C_result);	
	double C_rhov12;
	XL_readNumCell(XL_rhov12,C_rhov12,				" ARM_ERR: rhos12 numeric expected",C_result);	
	double C_rhov13;
	XL_readNumCell(XL_rhov13,C_rhov13,				" ARM_ERR: rhos13 numeric expected",C_result);	
	double C_rhov23;
	XL_readNumCell(XL_rhov23,C_rhov23,				" ARM_ERR: rhov23 numeric expected",C_result);	

	double C_rhoc12;
	XL_readNumCell(XL_rhoc12,C_rhoc12,				" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc21;
	XL_readNumCell(XL_rhoc21,C_rhoc21,				" ARM_ERR: rhoc21 numeric expected",C_result);	
	double C_rhoc13;
	XL_readNumCell(XL_rhoc13,C_rhoc13,				" ARM_ERR: rhoc13 numeric expected",C_result);	
	double C_rhoc31;
	XL_readNumCell(XL_rhoc31,C_rhoc31,				" ARM_ERR: rhoc31 numeric expected",C_result);	
	double C_rhoc23;
	XL_readNumCell(XL_rhoc23,C_rhoc23,				" ARM_ERR: rhoc12 numeric expected",C_result);	
	double C_rhoc32;
	XL_readNumCell(XL_rhoc32,C_rhoc32,				" ARM_ERR: rhoc32 numeric expected",C_result);	

	
	/// call the function
		long retCode=ARMLOCAL_Util_TriSABR_Eigenvalues(
			 C_rho1,	  C_rho2,	   C_rho3,
			 C_rhos12,  C_rhos23,  C_rhos13,
			 C_rhov12,  C_rhov23,  C_rhov13,
			 C_rhoc12,  C_rhoc13,
			 C_rhoc21,  C_rhoc23,
			 C_rhoc31,  C_rhoc32,
			C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = 6;
			int nbCols = 1;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Util_TriSABR_Eigenvalues" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/***********************************************

            CIR Parameter Model

	Inputs :
		Vector of model params Id

***********************************************/

LPXLOPER Local_CIR_ModelParamsCreate_Common(
	LPXLOPER XL_ParamsIdVec,
	bool PersistentInXL ){
	
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		vector<long> C_paramsIdVec( C_paramsIds.size() );
		for(int i = 0; i < C_paramsIds.size(); ++i )  	C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		exportFunc1Arg< vector<long>  >  ourFunc( C_paramsIdVec,  ARMLOCAL_CIR_ModelParamsCreate );
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}

	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EQHWVS_ParamModelCreate_Common" )
	return (LPXLOPER) &XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_CIR_ModelParamsCreate");
	bool PersistentInXL = true;
	return Local_CIR_ModelParamsCreate_Common(
        XL_ModelParamsId,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_PXL_CIR_ModelParamsCreate");
	bool PersistentInXL = false;
	return Local_CIR_ModelParamsCreate_Common(
        XL_ModelParamsId,
        PersistentInXL );
}
/***********************************************

            CIR Bond Price

	Inputs :
		return the bond price

***********************************************/

LPXLOPER Local_CIR_BondPrice_Common(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_phi,
	bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelParamId;
		XL_readStrCell(XL_ModelParamsId,C_ModelParamId," ARM_ERR: security id: object expected",C_result);

		double C_time;
		XL_readNumCell(XL_time,C_time," ARM_ERR: time : numeric expected",C_result);

		double D_r0=1.0;
		double C_r0;
		XL_readNumCellWD(XL_r0,C_r0,D_r0," ARM_ERR: r0 : numeric expected",C_result);

		double D_phi=1.0;
		double C_phi;
		XL_readNumCellWD(XL_phi,C_phi,D_phi," ARM_ERR: r0 : numeric expected",C_result);



		long retCode = ARMLOCAL_CIRBondPrice (LocalGetNumObjectId (C_ModelParamId),C_time, C_r0, C_phi, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}

	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondPrice(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_phi)
{
	ADD_LOG("Local_CIR_BondPrice");
	bool PersistentInXL = true;
	return Local_CIR_BondPrice_Common(
		XL_ModelParamsId,
		XL_time,
		XL_r0,
		XL_phi,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_BondPrice(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,	
	LPXLOPER XL_phi)
{
	ADD_LOG("Local_PXL_CIR_BondPrice");
	bool PersistentInXL = false;
	return Local_CIR_BondPrice_Common(
		XL_ModelParamsId,
		XL_time,
		XL_r0,
		XL_phi,
        PersistentInXL );
}
/***********************************************

            CIR Bond Density

	Inputs :
		return the bond density
***********************************************/

LPXLOPER Local_CIR_BondDensity_Common(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_r,
	LPXLOPER XL_frequency,
	bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelParamId;
		XL_readStrCell(XL_ModelParamsId,C_ModelParamId," ARM_ERR: security id: object expected",C_result);

		double C_time;
		XL_readNumCell(XL_time,C_time," ARM_ERR: time : numeric expected",C_result);

		double D_r0=1.0;
		double C_r0;
		XL_readNumCellWD(XL_r0,C_r0,D_r0," ARM_ERR: r0 : numeric expected",C_result);

		double C_r;
		XL_readNumCell(XL_r,C_r,," ARM_ERR: r : numeric expected",C_result);

		double D_frequency=1.0;
		double C_frequency;
		XL_readNumCellWD(XL_frequency,C_frequency,D_frequency," ARM_ERR: r0 : numeric expected",C_result);

		long retCode = ARMLOCAL_CIRBondDensity (LocalGetNumObjectId (C_ModelParamId),C_time, C_r0, C_r, C_frequency, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}

	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondDensity(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_r,
	LPXLOPER XL_frequency)
{
	ADD_LOG("Local_CIR_BondDensity");
	bool PersistentInXL = true;
	return Local_CIR_BondDensity_Common(
		XL_ModelParamsId,
		XL_time,
		XL_r0,
		XL_r,
		XL_frequency,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_BondDensity(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_r,
	LPXLOPER XL_frequency)
{
	ADD_LOG("Local_PXL_CIR_BondDensity");
	bool PersistentInXL = false;
	return Local_CIR_BondDensity_Common(
		XL_ModelParamsId,
		XL_time,
		XL_r0,
		XL_r,
		XL_frequency,
        PersistentInXL );
}
/***********************************************

            CIR Bond Distribution

	Inputs :
		return the bond distribution
***********************************************/

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondDistribution_Common(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_nbDiscr,
	LPXLOPER XL_factor,
	bool PersistentInXL ){
		ADD_LOG("Local_CIR_BondDistribution_Common");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	
		
		/// gets the inputs
		CCString C_ModelParamId;
		XL_readStrCell(XL_ModelParamsId,C_ModelParamId," ARM_ERR: security id: object expected",C_result);
	
		double C_time;
		XL_readNumCell(XL_time, C_time," ARM_ERR: time numeric expected",	C_result);
	
		double D_r0=1.0;
		double C_r0;
		XL_readNumCellWD(XL_r0,C_r0,D_r0," ARM_ERR: r0 : numeric expected",C_result);

		double C_nbDiscr;
		XL_readNumCell(XL_nbDiscr,C_nbDiscr," ARM_ERR: nbDiscr numeric expected",C_result);

		double D_factor=1.0;
		double C_factor;
		XL_readNumCellWD(XL_factor,C_factor,D_factor," ARM_ERR: factor : numeric expected",C_result);

		/// call the function
		long retCode=ARMLOCAL_CIRBondDistribution( LocalGetNumObjectId (C_ModelParamId), C_time, C_r0, C_nbDiscr, C_factor, C_result);
		
		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			int nbRows = (int) C_nbDiscr ;
			int nbCols = 3;
			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[nbCols*i+j] = C_result.getArray(nbCols*i+j);
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT
	
	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GaussianIntegrals_Legendre_Coeffs" )
	
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondDistribution(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_nbDiscr,
	LPXLOPER XL_factor)
{
	ADD_LOG("Local_CIR_BondDistribution");
	bool PersistentInXL = true;
	return Local_CIR_BondDistribution_Common(
		XL_ModelParamsId,
		XL_time,
		XL_r0,
		XL_nbDiscr,
		XL_factor,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_BondDistribution(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_nbDiscr,
	LPXLOPER XL_factor)
{
	ADD_LOG("Local_PXL_CIR_BondDistribution");
	bool PersistentInXL = false;
	return Local_CIR_BondDistribution_Common(
		XL_ModelParamsId,
		XL_time,
		XL_r0,
		XL_nbDiscr,
		XL_factor,
        PersistentInXL );
}

/***********************************************

            CIR Bond Weight

	Inputs :
		return the weight of the distribution
***********************************************/

LPXLOPER Local_CIR_Bond_Weight(
	LPXLOPER XL_ModelParamsId){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		bool PersistentInXL = true;
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelParamId;
		XL_readStrCell(XL_ModelParamsId,C_ModelParamId," ARM_ERR: security id: object expected",C_result);


		long retCode = ARMLOCAL_CIRBondWeight (LocalGetNumObjectId (C_ModelParamId),C_result);

		if(retCode == ARM_OK){
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}	
		else{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}

/***********************************************

            CIR Bond Expectation

	Inputs :
		return the weight of the distribution
***********************************************/

LPXLOPER Local_CIR_Bond_Expectation(
	LPXLOPER XL_ModelParamsId){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelParamId;
		XL_readStrCell(XL_ModelParamsId,C_ModelParamId," ARM_ERR: security id: object expected",C_result);

		long retCode = ARMLOCAL_CIRBondExpectation (LocalGetNumObjectId (C_ModelParamId),C_result);

		if(retCode == ARM_OK){
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}	
		else{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}

/***********************************************

            CIR Bond Variance

	Inputs :
		return the weight of the distribution
***********************************************/

LPXLOPER Local_CIR_Bond_Variance(
	LPXLOPER XL_ModelParamsId){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelParamId;
		XL_readStrCell(XL_ModelParamsId,C_ModelParamId," ARM_ERR: security id: object expected",C_result);

		long retCode = ARMLOCAL_CIRBondVariance (LocalGetNumObjectId (C_ModelParamId),C_result);

		if(retCode == ARM_OK){
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}	
		else{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}


/***********************************************

            Exp Riccati Create

	Inputs :
		build exponential riccati equation with y(t0)=y0

***********************************************/

__declspec(dllexport) LPXLOPER WINAPI Local_EXP_RICCATI_Create(
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_delta,
	LPXLOPER XL_lambda,
	LPXLOPER XL_x0,
	LPXLOPER XL_x1,
	LPXLOPER XL_x2,
	LPXLOPER XL_y0,
	LPXLOPER XL_t0 ){
	ADD_LOG("Local_EXP_RICCATI_Create");

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		double C_alpha;
		double C_beta;
		double C_delta;
		double C_lambda;
		double C_x0;
		double C_x1;
		double C_x2;
		double C_y0;
		double C_t0;

		XL_readNumCell(XL_alpha,	C_alpha,	" ARM_ERR: time : numeric expected",C_result);
		XL_readNumCell(XL_beta,		C_beta,		" ARM_ERR: time : numeric expected",C_result);
		XL_readNumCell(XL_delta,	C_delta,	" ARM_ERR: time : numeric expected",C_result);
		XL_readNumCell(XL_lambda,	C_lambda,	" ARM_ERR: time : numeric expected",C_result);

		XL_readNumCell(XL_x0,		C_x0,		" ARM_ERR: time : numeric expected",C_result);
		XL_readNumCell(XL_x1,		C_x1,		" ARM_ERR: time : numeric expected",C_result);
		XL_readNumCell(XL_x2,		C_x2,		" ARM_ERR: time : numeric expected",C_result);
	
		XL_readNumCell(XL_y0,		C_y0,		" ARM_ERR: time : numeric expected",C_result);	
		XL_readNumCell(XL_t0,		C_t0,		" ARM_ERR: time : numeric expected",C_result);	


		exportFunc9Args< double, double, double,
						 double, double, double,
						 double, double, double>
							ourFunc(	C_alpha, 	C_beta,		C_delta, 
										C_lambda, 	C_x0,		C_x1,
										C_x2, 		C_y0, 		C_t0,  ARMLOCAL_EXP_RICCATI_Create );

		fillXL_Result_withName(ourFunc, C_result, XL_result, true);
	}

	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EXP_RICCATI_Create" )
	return (LPXLOPER) &XL_result;
}


/***********************************************

            Exp Riccati Price

	Inputs :
		return the solution at time t of an exponential riccati equation

***********************************************/

__declspec(dllexport) LPXLOPER WINAPI Local_EXP_RICCATI_Price(
	LPXLOPER XL_RiccatiEq,
	LPXLOPER XL_t ){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_RiccatiEq;
		XL_readStrCell(XL_RiccatiEq, C_RiccatiEq," ARM_ERR: object expected",C_result);

		double C_t;
		XL_readNumCell(XL_t,C_t," ARM_ERR: numeric expected",C_result);	
	
		long retCode = ARMLOCAL_EXP_RICCATI_Price(LocalGetNumObjectId (C_RiccatiEq), C_t, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}

	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EXP_RICCATI_Int" )
	return (LPXLOPER)&XL_result;
}

/***********************************************

            Exp Riccati Int

	Inputs :
		return the integral of the solution of an exponential riccati equation between  t and T

***********************************************/

__declspec(dllexport) LPXLOPER WINAPI Local_EXP_RICCATI_Int(
	LPXLOPER XL_RiccatiEq,
	LPXLOPER XL_t,
	LPXLOPER XL_T ){

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_RiccatiEq;
		XL_readStrCell(XL_RiccatiEq, C_RiccatiEq," ARM_ERR: object expected",C_result);

		double C_t;
		double C_T;
		XL_readNumCell(XL_t, C_t," ARM_ERR: numeric expected",C_result);	
		XL_readNumCell(XL_T, C_T," ARM_ERR: numeric expected",C_result);	
		
		long retCode = ARMLOCAL_EXP_RICCATI_Int(LocalGetNumObjectId (C_RiccatiEq), C_t, C_T, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}

	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EXP_RICCATI_Int" )
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////////////////////////////
///
///	NonParametric Complete Option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_CompleteOption(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_Strike2_vec,
	LPXLOPER XL_Vol2_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S2params_vec,
	LPXLOPER XL_correlation,
	LPXLOPER XL_maturity,
	LPXLOPER XL_a1,
	LPXLOPER XL_b1,
	LPXLOPER XL_k1,
	LPXLOPER XL_a2,
	LPXLOPER XL_b2,
	LPXLOPER XL_k2,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algo,
	LPXLOPER XL_smiletype
	)
		{
			ADD_LOG("Local_Nonparametric_CompleteOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_Strike2_vec;
	XL_readNumVector(XL_Strike2_vec,C_Strike2_vec," ARM_ERR: initialparams Vector Strike2_vec: array of numeric expected",C_result);
	vector<double> C_Vol2_vec;
	XL_readNumVector(XL_Vol2_vec,C_Vol2_vec," ARM_ERR: initialparams Vector Vol2_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	vector<double> C_S2params_vec;
	XL_readNumVector(XL_S2params_vec,C_S2params_vec," ARM_ERR: initialparams Vector S2params_vec: array of numeric expected",C_result);
	double C_correlation;
	XL_readNumCell(XL_correlation,C_correlation,			" ARM_ERR: correlation numeric expected",C_result);	
	double C_maturity;
	XL_readNumCell(XL_maturity,C_maturity,			" ARM_ERR: maturity numeric expected",C_result);	
	double C_a1;
	XL_readNumCell(XL_a1,C_a1,			" ARM_ERR: a1 numeric expected",C_result);	
	double C_b1;
	XL_readNumCell(XL_b1,C_b1,			" ARM_ERR: b1 numeric expected",C_result);	
	double C_k1;
	XL_readNumCell(XL_k1,C_k1,			" ARM_ERR: k1 numeric expected",C_result);	

	double C_a2;
	XL_readNumCell(XL_a2,C_a2,			" ARM_ERR: a2 numeric expected",C_result);	
	double C_b2;
	XL_readNumCell(XL_b2,C_b2,			" ARM_ERR: b2 numeric expected",C_result);	
	double C_k2;
	XL_readNumCell(XL_k2,C_k2,			" ARM_ERR: k2 numeric expected",C_result);	
	double C_nbsteps;
	XL_readNumCell(XL_nbsteps,C_nbsteps,			" ARM_ERR: nbsteps numeric expected",C_result);	
	double C_algo;
	XL_readNumCell(XL_algo,C_algo,			" ARM_ERR: algo numeric expected",C_result);	
	double C_smiletype;
	XL_readNumCell(XL_smiletype,C_smiletype,	" ARM_ERR: smiletype numeric expected",C_result);	



	/// call the function
	long retCode=ARMLOCAL_Nonparametric_CompleteOption(
		C_Strike1_vec,
		C_Vol1_vec,
		C_Strike2_vec,
		C_Vol2_vec,
		C_S1params_vec,
		C_S2params_vec,
		C_correlation,
		C_maturity,
		C_a1,
		C_b1,
		C_k1,
		C_a2,
		C_b2,
		C_k2,
		C_nbsteps,
		C_algo,
		C_smiletype,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_CompleteOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_LogVolatility(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_k
	)
		{
			ADD_LOG("Local_Nonparametric_LogVolatility");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Nonparametric_LogVolatility(
		C_Strike1_vec,
		C_Vol1_vec,
		C_S1params_vec,
		C_k,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_LogVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_NormalVolatility(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_k
	)
		{
			ADD_LOG("Local_Nonparametric_NormalVolatility");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Nonparametric_NormalVolatility(
		C_Strike1_vec,
		C_Vol1_vec,
		C_S1params_vec,
		C_k,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_NormalVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_NormalDistribution(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	)
		{
			ADD_LOG("Local_Nonparametric_NormalDistribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	
	double C_S;
	XL_readNumCell(XL_S,C_S,			" ARM_ERR: S numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Nonparametric_NormalDistribution(
		C_Strike1_vec,
		C_Vol1_vec,
		C_S1params_vec,
		C_S,
		C_T,
		C_k,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_NormalDistribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_NormalQuantile(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	)
		{
			ADD_LOG("Local_Nonparametric_NormalQuantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	
	double C_S;
	XL_readNumCell(XL_S,C_S,			" ARM_ERR: S numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Nonparametric_NormalQuantile(
		C_Strike1_vec,
		C_Vol1_vec,
		C_S1params_vec,
		C_S,
		C_T,
		C_k,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_NormalQuantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_LogNormalDistribution(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	)
		{
			ADD_LOG("Local_Nonparametric_LogNormalDistribution");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	
	double C_S;
	XL_readNumCell(XL_S,C_S,			" ARM_ERR: S numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Nonparametric_LogNormalDistribution(
		C_Strike1_vec,
		C_Vol1_vec,
		C_S1params_vec,
		C_S,
		C_T,
		C_k,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_LogNormalDistribution" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_LogNormalQuantile(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	)
		{
			ADD_LOG("Local_Nonparametric_LogNormalQuantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	vector<double> C_Strike1_vec;
	XL_readNumVector(XL_Strike1_vec,C_Strike1_vec," ARM_ERR: initialparams Vector Strike1_vec: array of numeric expected",C_result);
	vector<double> C_Vol1_vec;
	XL_readNumVector(XL_Vol1_vec,C_Vol1_vec," ARM_ERR: initialparams Vector Vol1_vec: array of numeric expected",C_result);
	vector<double> C_S1params_vec;
	XL_readNumVector(XL_S1params_vec,C_S1params_vec," ARM_ERR: initialparams Vector S1params_vec: array of numeric expected",C_result);
	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	
	double C_S;
	XL_readNumCell(XL_S,C_S,			" ARM_ERR: S numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	


	/// call the function
	long retCode=ARMLOCAL_Nonparametric_LogNormalQuantile(
		C_Strike1_vec,
		C_Vol1_vec,
		C_S1params_vec,
		C_S,
		C_T,
		C_k,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Nonparametric_LogNormalQuantile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BiShiftedHeston_VanillaOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_V1,
	LPXLOPER XL_Vinfini1,
	LPXLOPER XL_lambda1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_gamma1,
	LPXLOPER XL_F2,
	LPXLOPER XL_V2,
	LPXLOPER XL_Vinfini2,
	LPXLOPER XL_lambda2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_gamma2,
	LPXLOPER XL_Correlations_vec,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_LambdaB,
	LPXLOPER XL_Flag
	)
		{
			ADD_LOG("Local_BiShiftedHeston_VanillaOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

	/// gets the inputs
	double C_F1;
	XL_readNumCell(XL_F1,C_F1,			" ARM_ERR: F1 numeric expected",C_result);	
	double C_V1;
	XL_readNumCell(XL_V1,C_V1,			" ARM_ERR: V1 numeric expected",C_result);	
	double C_Vinfini1;
	XL_readNumCell(XL_Vinfini1,C_Vinfini1,			" ARM_ERR: Vinfini1 numeric expected",C_result);	
	double C_lambda1;
	XL_readNumCell(XL_lambda1,C_lambda1,			" ARM_ERR: lambda1 numeric expected",C_result);	
	double C_nu1;
	XL_readNumCell(XL_nu1,C_nu1,			" ARM_ERR: nu1 numeric expected",C_result);	
	double C_rho1;
	XL_readNumCell(XL_rho1,C_rho1,			" ARM_ERR: rho1 numeric expected",C_result);	
	double C_gamma1;
	XL_readNumCell(XL_gamma1,C_gamma1,			" ARM_ERR: gamma1 numeric expected",C_result);	

	double C_F2;
	XL_readNumCell(XL_F2,C_F2,			" ARM_ERR: F2 numeric expected",C_result);	
	double C_V2;
	XL_readNumCell(XL_V2,C_V2,			" ARM_ERR: V2 numeric expected",C_result);	
	double C_Vinfini2;
	XL_readNumCell(XL_Vinfini2,C_Vinfini2,			" ARM_ERR: Vinfini2 numeric expected",C_result);	
	double C_lambda2;
	XL_readNumCell(XL_lambda2,C_lambda2,			" ARM_ERR: lambda2 numeric expected",C_result);	
	double C_nu2;
	XL_readNumCell(XL_nu2,C_nu2,			" ARM_ERR: nu2 numeric expected",C_result);	
	double C_rho2;
	XL_readNumCell(XL_rho2,C_rho2,			" ARM_ERR: rho2 numeric expected",C_result);	
	double C_gamma2;
	XL_readNumCell(XL_gamma2,C_gamma2,			" ARM_ERR: gamma2 numeric expected",C_result);	

	vector<double> C_Correlations_vec;
	XL_readNumVector(XL_Correlations_vec,C_Correlations_vec," ARM_ERR: Correlations Vector : array of numeric expected",C_result);

	double C_k;
	XL_readNumCell(XL_k,C_k,			" ARM_ERR: k numeric expected",C_result);	
	double C_T;
	XL_readNumCell(XL_T,C_T,			" ARM_ERR: T numeric expected",C_result);	
	double C_CallPut;
	XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: call or put string expected",C_result);
	double C_LambdaB;
	XL_readNumCell(XL_LambdaB,C_LambdaB,			" ARM_ERR: LambdaB numeric expected",C_result);	
	double C_Flag;
	XL_readNumCell(XL_Flag,C_Flag,			" ARM_ERR: Flag numeric expected",C_result);	

	/// call the function
	long retCode=ARMLOCAL_BiShiftedHeston_VanillaOption(
		C_F1,
		C_V1,
		C_Vinfini1,
		C_lambda1,
		C_nu1,
		C_rho1,
		C_gamma1,
		C_F2,
		C_V2,
		C_Vinfini2,
		C_lambda2,
		C_nu2,
		C_rho2,
		C_gamma2,
		C_Correlations_vec,
		C_k,
		C_T,
		C_CallPut,
		C_LambdaB,
		C_Flag,
		C_result);

	/// return the result as an LPXLOPER
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
	}
	else
	{
		ARM_ERR();
	}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiShiftedHeston_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport)LPXLOPER WINAPI Local_Merton_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Sigma,
	LPXLOPER XL_Lambda1,
	LPXLOPER XL_U1,
	LPXLOPER XL_Lambda2,
	LPXLOPER XL_U2,
	LPXLOPER XL_N)
{
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN

	{/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		double  fwd, k, t, callput, sigma, lambda1, U1, lambda2, U2, N, defN = 20.;

		XL_readNumCell(XL_F, fwd, "ARM_ERR : fwd: numeric expeted", C_result);
		XL_readNumCell(XL_K, k, "ARM_ERR : strike : numeric expeted", C_result);
		XL_readNumCell(XL_T, t, "ARM_ERR : maturity : numeric expeted", C_result);
		XL_readNumCell(XL_Sigma, sigma, "ARM_ERR : sigma: numeric expeted", C_result);
		XL_readNumCell(XL_Lambda1, lambda1, "ARM_ERR : lambda1: numeric expeted", C_result);
		XL_readNumCell(XL_Lambda2, lambda2, "ARM_ERR : lambda2: numeric expeted", C_result);
		XL_readNumCell(XL_U1, U1, "ARM_ERR : U1: numeric expeted", C_result);
		XL_readNumCell(XL_U2, U2, "ARM_ERR : U2: numeric expeted", C_result);
		XL_readNumCellWD(XL_N, N, defN, "ARM_ERR : N : numeric expected", C_result);
		XL_GETCONVCALLORPUT(XL_CallPut,callput," ARM_ERR: cap or floor string expected",C_result);

		long retCode = ARMLOCAL_MertonOption(fwd, k, t, callput, sigma, lambda1, U1, lambda2, U2, N, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Merton_VanillaOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport)LPXLOPER WINAPI Local_BSImpliedVol(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Target,
	LPXLOPER XL_BondPrice)
{
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN

	{/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		double fwd, k, t, callput, target, bp, def(1.);

		XL_readNumCell(XL_F, fwd, "ARM_ERR : fwd: numeric expeted", C_result);
		XL_readNumCell(XL_K, k, "ARM_ERR : strike : numeric expeted", C_result);
		XL_readNumCell(XL_T, t, "ARM_ERR : maturity : numeric expeted", C_result);
		XL_readNumCell(XL_Target, target, "ARM_ERR : target: numeric expeted", C_result);
		XL_readNumCellWD(XL_BondPrice, bp, def, "ARM_ERR : bond price : numeric expected", C_result);
		XL_GETCONVCALLORPUT(XL_CallPut,callput," ARM_ERR: cap or floor string expected",C_result);

		long retCode = ARMLOCAL_BSImpliedVol(fwd, k, t, callput, target/bp, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSImpliedVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport)LPXLOPER WINAPI Local_SuperNormal_Heston_VanillaOption(
	LPXLOPER XL_S0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Theta1,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_Nu1,
	LPXLOPER XL_V01,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Theta2,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_V02,
	LPXLOPER XL_Nb)
{
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN

	{/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		double C_S0,C_K,C_T,C_CallPut,C_V01,C_Rho1,C_Theta1,C_Kappa1,C_Nu1,C_V02,C_Rho2,C_Theta2,C_Kappa2,C_Nu2,C_nb, defNb(20.);

		XL_readNumCell(XL_S0, C_S0, "ARM_ERR : fwd: numeric expeted", C_result);
		XL_readNumCell(XL_K, C_K, "ARM_ERR : strike : numeric expeted", C_result);
		XL_readNumCell(XL_T, C_T, "ARM_ERR : maturity : numeric expeted", C_result);
		XL_GETCONVCALLORPUT(XL_CallPut,C_CallPut," ARM_ERR: cap or floor string expected",C_result);
		XL_readNumCell(XL_Rho1, C_Rho1, "ARM_ERR : rho1: numeric expeted", C_result);
		XL_readNumCell(XL_Kappa1, C_Kappa1, "ARM_ERR : kappa1: numeric expeted", C_result);
		XL_readNumCell(XL_Theta1, C_Theta1, "ARM_ERR : theta1: numeric expeted", C_result);
		XL_readNumCell(XL_Nu1, C_Nu1, "ARM_ERR : nu1: numeric expeted", C_result);
		XL_readNumCell(XL_V01, C_V01, "ARM_ERR : V01: numeric expeted", C_result);
		XL_readNumCell(XL_Rho2, C_Rho2, "ARM_ERR : rho2: numeric expeted", C_result);
		XL_readNumCell(XL_Kappa2, C_Kappa2, "ARM_ERR : kappa2: numeric expeted", C_result);
		XL_readNumCell(XL_Theta2, C_Theta2, "ARM_ERR : theta2: numeric expeted", C_result);
		XL_readNumCell(XL_Nu2, C_Nu2, "ARM_ERR : nu2: numeric expeted", C_result);
		XL_readNumCell(XL_V02, C_V02, "ARM_ERR : V02: numeric expeted", C_result);
		XL_readNumCellWD(XL_Nb, C_nb, defNb, "ARM_ERR : nb : numeric expected", C_result);

		long retCode = ARMLOCAL_SuperNormal_Heston_VanillaCall(C_Rho1,C_Kappa1,C_Theta1,C_Nu1,C_V01,
															   C_Rho2,C_Kappa2,C_Theta2,C_Nu2,C_V02,
															   C_S0,C_K,C_T,C_CallPut,C_nb, C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSImpliedVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
