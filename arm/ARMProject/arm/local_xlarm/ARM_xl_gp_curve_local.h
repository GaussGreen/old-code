/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_curves_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARM_XL_GP_CURVE_LOCAL_H
#define ARM_XL_GP_CURVE_LOCAL_H


#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

//////////////////////////////////
/// 
/// In order to factorise code
/// many function are calling
/// common function
/// these functions are not
/// declared here 
/// 
/// only exported functions are 
/// included here to facilitate 
/// the reading of the header
/// 
//////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_Create(
	LPXLOPER XL_abscisses,
	LPXLOPER XL_ordinates,
	LPXLOPER XL_interpolatorType,
	LPXLOPER XL_sortAbscisses,
	LPXLOPER XL_alwaysMulti);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenericCurve_Create(
	LPXLOPER XL_abscisses,
	LPXLOPER XL_ordinates,
	LPXLOPER XL_interpolatorType,
	LPXLOPER XL_sortAbscisses,
	LPXLOPER XL_alwaysMulti);

__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_Interpolate(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_abscisse );

__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_CptCurve(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_NewAbscisses );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenericCurve_CptCurve(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_NewAbscisses );

__declspec(dllexport) LPXLOPER WINAPI Local_GenericCurve_Insert(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_NewAbscisse,
	LPXLOPER XL_NewOrdinate );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenericCurve_Insert(
	LPXLOPER XL_CurveId,
	LPXLOPER XL_NewAbscisse,
	LPXLOPER XL_NewOrdinate );

///////////////////////////////////
/// addin to create a flat surface
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FlatSurface_Create(
	LPXLOPER XL_Value);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FlatSurface_Create(
	LPXLOPER XL_Value);

///////////////////////////////////
/// addin to create a linear surface
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LinSurface_Create(
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3,
    LPXLOPER XL_InterpoleType);

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LinSurface_Create(
	LPXLOPER XL_X1,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3,
    LPXLOPER XL_InterpoleType);

///////////////////////////////////
/// addin to interpolate a surface
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Surface_Interpolate(
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_X2,
	LPXLOPER XL_X3);


///////////////////////////////////
/// addin to convert vol summit to linear surface
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FromVolSummitToSurfce_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_InterpolType);

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FromVolSummitToSurfce_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_InterpolType);

///////////////////////////////////
/// addin to convert vol summit to ARM_Curve
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FromVolSummitToCurve_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_calcMethod);

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FromVolSummitToCurve_Create(
	LPXLOPER XL_VolId,
    LPXLOPER XL_calcMethod);

#endif