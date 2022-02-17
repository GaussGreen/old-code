/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_curve.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_CURVE_H
#define ARMLOCAL_GP_CURVE_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_result.h"
#include <vector>
using std::vector;


////////////////////////////////////////////
//// function to create a generic curve
////////////////////////////////////////////
extern long ARMLOCAL_GenericCurve_Create(
	const vector<double>&	C_abscisses,
	const vector<double>&   C_ordinates,
	const long&				C_rowsNb,
	const long&				C_colsNb,
	const bool&				C_sortAbscisses,
	const string&			C_interpolatorName,
	const bool&				C_alwaysMulti,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// function to interpolate a generic curve
////////////////////////////////////////////
extern long ARMLOCAL_GenericCurve_Interpolate(
	long					C_genericCurveId,
	double					ordinate,
	vector<double>&			vecResult,
    ARM_result&				result );


////////////////////////////////////////////
//// function to cpt a curve from another generic curve
////////////////////////////////////////////
extern long ARMLOCAL_GenericCurve_CptCurve(
	long					C_genericCurveId,
	const vector<double>&	C_newAbscisses,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// function to cpt a curve from another generic curve
////////////////////////////////////////////
extern long ARMLOCAL_GenericCurve_Insert(
	long					C_genericCurveId,
	double					C_newAbscisse,
	const vector<double>&   C_newOrdinate,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID );



////////////////////////////////////////////
//// function to create a flat curve
////////////////////////////////////////////

extern long ARMLOCAL_FlatSurface_Create(
	const double&	value,
	ARM_result&		result, 
    long			objId = ARM_NULL_OBJECT_ID );


////////////////////////////////////////////
//// function to create a linear 2D surface
////////////////////////////////////////////
extern long ARMLOCAL_LinSurface_Create(
	const vector<double>& X1,
	const vector<double>& X2,
	const vector<double>& X3,
	const long& X3RowNb,
	const long& X3ColNb,
    const string& type,
	ARM_result&	result, 
	long		objId = ARM_NULL_OBJECT_ID );


////////////////////////////////////////////
//// function to interpolate a surface
////////////////////////////////////////////
extern long ARMLOCAL_Surface_Interpolate(
	const long& SurfaceId,
	const double& X1,
	const double& X2,
	ARM_result&	result, 
	long		objId = ARM_NULL_OBJECT_ID );


////////////////////////////////////////////
//// function to insert a point in a surface
////////////////////////////////////////////
extern long ARMLOCAL_Surface_Insert(
	const long& surfaceId,
	const double& x1,
	const double& x2,
	const double& x3,
	ARM_result&	result, 
	long		objId = ARM_NULL_OBJECT_ID );


////////////////////////////////////////////
//// function to convert summit volatility to linear surface
////////////////////////////////////////////
extern long ARMLOCAL_FromVolSummitToSurfce_Create(
	const long& volId,
    const string& typeStr,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// function to create a correl curve matrix
////////////////////////////////////////////
extern long ARMLOCAL_CurveMatrix_Create(
	const VECTOR<long>& correlCurvesId,
    const long& nbRows,
	const long& nbCols,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// function to convert summit volatility to ARM_Curve
////////////////////////////////////////////
extern long ARMLOCAL_FromVolSummitToCurve_Create(
	const long& volId,
    const long& calcMethod,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID );

#endif
