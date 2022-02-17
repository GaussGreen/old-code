/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculators.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */


#include "firstToBeIncluded.h"
#include "ARM_local_gp_curve.h"

#include "ARM_local_wrapper.h"
#include "ARM_local_glob.h"
#include "ARM_local_class.h"

/// gpbase
#include <GP_Base\gpbase\curvetypedef.h>
#include <GP_Base\gpbase\curvefactory.h>
#include <GP_Base\gpbase\curve.h>
#include <GP_Base\gpbase\curvematrix.h>
#include <GP_Base\gpbase\vectorarithmetic.h>
#include <GP_Base\gpbase\vectormanip.h>
#include <GP_Base\gpbase\gpvector.h>
#include <GP_Base\gpbase\gpmatrix.h>
#include <GP_Base\gpbase\surface.h>
#include <GP_Base\gpbase\interpolator2D.h>
#include <GP_Base\gpbase\surfacetypedef.h>
#include <GP_Base\gpbase\gplinalgtypedef.h>
#include <GP_Base\gpbase\typedef.h>
#include <GP_Base\gpbase\argconvdefault.h>
#include <GP_Base\gpbase\gplinalgconvert.h>

/// gphelp
#include <GP_Help\gphelp\crmcookies.h>


/// Kernel
#include "volcurv.h"

using ARM::ARM_Curve;
using ARM::ARM_LinInterpCstExtrapolDble;
using ARM::ARM_StepUpRightOpenCstExtrapolDble;
using ARM::ARM_StepUpLeftOpenCstExtrapolDble;
using ARM::ARM_MultiCurve;
using ARM::ARM_LinInterpCstExtrapolVec;
using ARM::ARM_StepUpRightOpenCstExtrapolVec;
using ARM::ARM_StepUpLeftOpenCstExtrapolVec;
using ARM::ARM_FlatSurface;
using ARM::ARM_SurfaceWithInterpol;
using ARM::ARM_Interporlator2D;
using ARM::ARM_2DLin2Interpol;
using ARM::ARM_2DLin1Interpol;
using ARM::ARM_GP_T_Vector;
using ARM::ARM_GP_Vector;
using ARM::ARM_GP_T_Matrix;
using ARM::ARM_GP_Matrix;
using ARM::To_ARM_GP_Matrix;
using ARM::To_ARM_GP_Vector;
using ARM::ARM_Surface;
using ARM::ARM_CurveFactory;
using ARM::ARM_GP_Vector;
using ARM::ARM_ArgConv_InterpolType;
using ARM::ARM_ArgConvReverse_InterpolType;
using ARM::ARM_InterpolationType; 
using ARM::ARM_InterpolType; 
using ARM::ARM_CRMCookies;
using ARM::ARM_CurveMatrix;
using ARM::ARM_Interpolator;

////////////////////////////////////////////
//// Function to create a genericCurve
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
    long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* genericCurve = NULL;

	try
	{

		genericCurve = ARM_CurveFactory::CreateGenericCurve(
			C_abscisses,
			C_ordinates,
			C_rowsNb,
			C_colsNb,
			C_sortAbscisses,
			C_interpolatorName,
			C_alwaysMulti);

		/// assign object
		if( !assignObject( genericCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete genericCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_GenericCurve_Interpolate(
	long			C_genericCurveId,
	double			abscisse,
	vector<double>&	vecResult,
    ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_genericCurveId);
		if( !LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( armObj, ARM_GENERIC_CURVE) )
		{
			result.setMsg ("ARM_ERR: generic genericCurve is not of a good type");
			return ARM_KO;
		};

		/// test type
		if( ARM_Curve*  simpleCurve = dynamic_cast< ARM_Curve* >( armObj ) )
		{
			vecResult = vector<double>( 1, simpleCurve->Interpolate( abscisse ) );
		}
		else if( ARM_MultiCurve* multiCurve = dynamic_cast< ARM_MultiCurve* >( armObj ) )
		{
			ARM_GP_Vector res = multiCurve->Interpolate( abscisse );
			vecResult = vector<double>( res.size() );
			std::copy(res.begin(), res.end(), vecResult.begin() );
		}
		else
		{
			CCString local_msg ("unknown curve: permitted are simple curve and multicurve");
			result.setMsg (local_msg);
			return ARM_KO;
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}




extern long ARMLOCAL_GenericCurve_CptCurve(
	long					C_genericCurveId,
	const vector<double>&	C_newAbscisses,
	ARM_result&				result, 
    long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* genericCurve = NULL;

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_genericCurveId);
		if( !LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( armObj, ARM_GENERIC_CURVE) )
		{
			result.setMsg ("ARM_ERR: generic genericCurve is not of a good type");
			return ARM_KO;
		};

		/// test type
		if( ARM_Curve*  simpleCurve = dynamic_cast< ARM_Curve* >( armObj ) )
		{
			genericCurve = simpleCurve->CptCurve( ARM_GP_Vector(C_newAbscisses));
		}
		else if( ARM_MultiCurve* multiCurve = dynamic_cast< ARM_MultiCurve* >( armObj ) )
		{
			genericCurve = multiCurve->CptCurve( ARM_GP_Vector(C_newAbscisses));
		}
		else
		{
			CCString local_msg ("unknown curve: permitted are simple curve and multicurve");
			result.setMsg (local_msg);
			return ARM_KO;
		}

		/// assign object
		if( !assignObject( genericCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete genericCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_OK;
}



extern long ARMLOCAL_GenericCurve_Insert(
	long					C_genericCurveId,
	double					C_newAbscisse,
	const vector<double>&   C_newOrdinate,
	ARM_result&				result, 
    long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* genericCurve = NULL;

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_genericCurveId);
		if( !LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( armObj, ARM_GENERIC_CURVE) )
		{
			result.setMsg ("ARM_ERR: generic genericCurve is not of a good type");
			return ARM_KO;
		};

		/// simple form
		if( ARM_Curve*  simpleCurve = dynamic_cast< ARM_Curve* >( armObj ) )
		{
			if( C_newOrdinate.size() != 1 )
			{
				result.setMsg ("ARM_ERR: a curve is expecting only a number for an ordinate!");
				return ARM_KO;
			};
			simpleCurve = (ARM_Curve*) simpleCurve->Clone();
			simpleCurve->insert(C_newAbscisse,C_newOrdinate[0] );
			genericCurve = simpleCurve;
		}

		/// vectorial form
		else if( ARM_MultiCurve* multiCurve = dynamic_cast< ARM_MultiCurve* >( armObj ) )
		{
			if( C_newOrdinate.size() != multiCurve[0].size() )
			{
				result.setMsg ("ARM_ERR: ordinate of wrong size!");
				return ARM_KO;
			};
			multiCurve = (ARM_MultiCurve*) multiCurve->Clone();
			multiCurve->insert(C_newAbscisse,ARM_GP_Vector(C_newOrdinate));
			genericCurve = multiCurve;
		}
		else
		{
			CCString local_msg ("unknown curve: permitted are simple curve and multicurve");
			result.setMsg (local_msg);
			return ARM_KO;
		}

		/// assign object
		if( !assignObject( genericCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete genericCurve;
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_OK;
}


////////////////////////////////////////////
//// Function to create a flat surface
////////////////////////////////////////////
extern long ARMLOCAL_FlatSurface_Create(
	const double& value,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_FlatSurface* surface= NULL;

	try
	{
		surface = new ARM_FlatSurface( value );

		/// assign object
		if( !assignObject( surface, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete surface;
		x.DebugPrint();
		ARM_RESULT();
	}
}




////////////////////////////////////////////
//// Function to create a linear surface
////////////////////////////////////////////
extern long ARMLOCAL_LinSurface_Create(
	const vector<double>& X1,
	const vector<double>& X2,
	const vector<double>& X3,
	const long& X3RowNb,
	const long& X3ColNb,
    const string& typeStr,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_SurfaceWithInterpol* surface		= NULL;

	if( X1.size() != X3RowNb )
	{
		CCString local_msg ("X1 size != X3 Rows");
		result.setMsg (local_msg);
		return ARM_KO;
	}

	if( X2.size() != X3ColNb )
	{
		CCString local_msg ("X2 size != X3 Cols");
		result.setMsg (local_msg);
		return ARM_KO;
	}
	
	try
	{
		ARM_GP_Vector x1Vec(X1);
		ARM_GP_Vector x2Vec(X2);
		ARM_GP_Matrix x3Mat(X3RowNb,X3ColNb,X3);
		ARM_InterpolType type = (ARM_InterpolType) ARM_ArgConv_InterpolType.GetNumber(typeStr);

		surface		= new ARM_SurfaceWithInterpol( x1Vec, x2Vec, x3Mat,type);

		/// assign object
		if( !assignObject( surface , result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete surface;
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// function to interpolate a surface
////////////////////////////////////////////
extern long ARMLOCAL_Surface_Interpolate(
	const long& C_SurfaceId,
	const double& X1,
	const double& X2,
	ARM_result&	result, 
	long		objId  )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(C_SurfaceId);
		ARM_Surface* surface	= dynamic_cast<ARM_Surface*>(armObj);

		if( !surface )
		{
			result.setMsg ("ARM_ERR: surface is not of a good type");
			return ARM_KO;
		}
		else
		{
			result.setDouble( surface->Interpolate(X1,X2) );
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a linear surface
////////////////////////////////////////////
extern long ARMLOCAL_Surface_Insert(
	const long& surfaceId,
	const double& x1,
	const double& x2,
	const double& x3,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Surface* oldSurface		= NULL;
	ARM_Surface* newSurface		= NULL;

	try
	{
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(surfaceId);
		oldSurface				= dynamic_cast<ARM_Surface*>(armObj);

		if( !oldSurface )
		{
			result.setMsg ("ARM_ERR: surface is not of a good type");
			return ARM_KO;
		}
		else
		{
			newSurface = (ARM_Surface*) oldSurface->Clone();
			newSurface->insert(x1,x2,x3);
		}	

		/// assign object
		if( !assignObject( newSurface, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete oldSurface;
		delete newSurface;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to convert vol summit to linear surface
////////////////////////////////////////////
extern long ARMLOCAL_FromVolSummitToSurfce_Create(
	const long& volId,
    const string& typeStr,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_VolCurve* volatility				= NULL;
	ARM_SurfaceWithInterpol* surface		= NULL;

	try
	{
			/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Summit vol converting to ARM_GP_Surface" );

		if( !GetObjectFromIdWithDynamicCastCheck( &volatility, volId  ) )
		{
			result.setMsg ("ARM_ERR: Volatility from summit is not of a good type");
			return ARM_KO;
		}
		ARM_Matrix* volMatrix = volatility->GetVolatilities();
		ARM_GP_Matrix x3Mat(To_ARM_GP_Matrix(*volMatrix));
		ARM_Vector* expiries = volatility->GetExpiryTerms();

		ARM_GP_Vector x1Vec(To_ARM_GP_Vector(*expiries));
		ARM_GP_Vector x2Vec(To_ARM_GP_Vector(*(volatility->GetStrikes())));

		ARM_InterpolType type = (ARM_InterpolType) ARM_ArgConv_InterpolType.GetNumber(typeStr);

		x1Vec*=365;
		x3Mat/=100.0;
		surface		= new ARM_SurfaceWithInterpol( x1Vec, x2Vec, x3Mat,type);

		/// assign object
		if( !assignObject( surface , result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete surface;
		x.DebugPrint();
		ARM_RESULT();
	}
}
////////////////////////////////////////////
//// Function to convert vol summit to linear surface
////////////////////////////////////////////
extern long ARMLOCAL_CurveMatrix_Create(
	const VECTOR<long>& correlCurvesId,
    const long& nbRows,
	const long& nbCols,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CurveMatrix* curveMatrix		= NULL;
	ARM_Curve* curve;
	vector<ARM_Curve*> curves;

	try
	{
			/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Creating Curve Correl Matrix" );

		for (size_t i = 0; i < correlCurvesId.size(); ++i)
		{
			if( !GetObjectFromIdWithDynamicCastCheck( &curve, correlCurvesId[i]  ) )
			{
				result.setMsg ("ARM_ERR: Correl Curve is not of a good type");
				return ARM_KO;
			}

			curves.push_back(curve);
		}
		
		curveMatrix = new ARM_CurveMatrix(curves,nbRows,nbCols);

		/// assign object
		if( !assignObject( curveMatrix , result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete curveMatrix;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to convert vol summit to ARM_Curve
////////////////////////////////////////////
extern long ARMLOCAL_FromVolSummitToCurve_Create(
	const long& volId,
    const long& calcMethod,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_VolCurve* volatility	= NULL;
	ARM_Curve* curve			= NULL;

	ARM_Interpolator<double,double>* interpolator = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Summit vol converting to ARM_Curve" );

		if( !GetObjectFromIdWithDynamicCastCheck( &volatility, volId  ) )
		{
			result.setMsg ("ARM_ERR: Volatility from summit is not of a good type");
			return ARM_KO;
		}

		ARM_GP_Vector abscisses = To_ARM_GP_Vector(volatility->GetExpiryTerms());
		ARM_GP_Vector ordinates = To_ARM_GP_Vector(volatility->GetVolatilities()->GetColumn(0));

		abscisses *= 365;

		switch (calcMethod)
		{
			case K_STEPUP_RIGHT:
					interpolator = new ARM::ARM_StepUpRightOpenCstExtrapolDble();
				break;
			case K_STEPUP_LEFT:
					interpolator = new ARM::ARM_StepUpLeftOpenCstExtrapolDble();
				break;
			default:
				interpolator = new ARM::ARM_LinInterpCstExtrapolDble();
		}


		curve		= new ARM_Curve( abscisses, ordinates, interpolator);

		/// assign object
		if( !assignObject( curve , result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete curve;
		x.DebugPrint();
		ARM_RESULT();
	}
}