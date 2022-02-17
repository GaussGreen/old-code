/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_util.cpp,v $
 * Revision 1.1  2004/01/13 15:08:43  jmprie
 * Initial version
 *
 */


/*! \file ARM_local_gp_util.cpp
 *
 *  \brief file for vairous utility function for the generic pricer local addins functions
 *
 *	\author  J.M. Prie
 *	\version 1.0
 *	\date January 2004
 */

#include "firstToBeIncluded.h"

///////////////////////////////////////////////
/// WARNING include these headers FIRST
/// because it uses the new style headers
/// CCxl uses cstyle headers and leads to conflict
/// if defined first
///////////////////////////////////////////////

#include <gpbase\curve.h>
#include <gpbase\surface.h>
#include <gpbase\surfacetypedef.h>
#include <gpbase\autocleaner.h>
#include <gpbase\eventviewerrep.h>
#include <gpbase\errviewer.h>
#include <gpbase\eventviewer.h>
#include <gpbase\gpmatrix.h>
#include <gpbase\gptrigomatrix.h>
#include <gpbase\gpvector.h>
#include <gpbase\interpolator2D.h>
#include <gpbase\gpmatrixlinalg.h>
#include <gpbase\typedef.h>
#include <gpbase\argconvdefault.h>

#include <gpinfra\typedef.h>
#include <gpinfra\modelparams.h>
#include <gpinfra\pricingmodel.h>
#include <gpinfra\pricingmodelir.h>
#include <GP_Help\gphelp\crmcookies.h>
#include <gpmodels\MultiAssets.h>
#include <gpmodels\modelparamshw.h>

#include <gpclosedforms\fxconvertdata.h>

#include <gpnumlib\typedef.h>
#include <gpnumlib\argconvdefault.h>
#include <gpnumlib\regression.h>

#include <gpclosedforms\argconvdefault_CF.h>

#include "ARM_result.h"
#include "ARM_local_glob.h"
#include "ARM_local_wrapper.h"

#include "glob\expt.h"

#include <gpinfra\gramfunctorarg.h>
#include <gpinfra\gramfunctorconv.h>

#include <mod\bsmodel.h>
#include <inst\portfolio.h>



//// using the namespace directive to access ARM object!
using ARM::ComputeEndDate;
using ARM::ARM_GramFctorArg;
using ARM::ARM_ModelParams;
using ARM::ARM_ModelParam;
using ARM::ARM_Curve;
using ARM::ARM_PricingModel;
using ARM::ARM_CRMCookies;
using ARM::ARM_PricingModelIR;
using ARM::ARM_ModelParamsHW;
using ARM::ARM_EventViewerRep;
using ARM::ARM_ErrViewer;
using ARM::ARM_EventViewerImp;
using ARM::ARM_GP_Matrix;
using ARM::ARM_GP_MatrixPtr;
using ARM::TrigoMatrix;
using ARM::std::vector<double>;
using ARM::ARM_GP_StrVector;
using ARM::ARM_GP_VectorPtr;
using ARM::ARM_AutoCleaner;
using ARM::LeastSquareRegression;
using ARM::ACPTransformation;
using ARM::ARM_Surface;
using ARM::ARM_SurfaceWithInterpol;
using ARM::ARM_2DLin1Interpol;
using ARM::ARM_LinInterpCstExtrapolDble;
using ARM::ARM_ZeroCurvePtr;
using ARM::ARM_MultiAssetsModel;
using ARM::ARM_ArgConv_InterpolType;
using ARM::ARM_ArgConvReverse_InterpolType;
using ARM::ARM_InterpolationType; 
using ARM::ARM_InterpolType; 
using ARM::ARM_Regression;
using ARM::ARM_RegressionPtr;
using ARM::ARM_LSRegression;
using ARM::ARM_LOESSRegression;
using ARM::ARM_ArgConv_RegMode;
using ARM::ARM_FXMktDataToTotemFormat;
using ARM::ARM_ArgConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag;


static const double varCovarStepTime = 0.0001;
#define XL_TYPE_STRING  2

////////////////////////////////////////////
//// XL Dates to GP Times convertion
////////////////////////////////////////////
long ConvertDatesToTimes(
   ARM_PricingModel& mod,
   double       fromDate,
   double       toDate,
   ARM_Date     startDate,
   ARM_Date     endDate,
   double&      fromTime,
   double&      toTime,
   double&      startTime,
   double&      endTime)
{
    /// Convert dates to times
    char fromDateStr[11];
	Local_XLDATE2ARMDATE(fromDate,fromDateStr);
    fromTime=mod.GetTimeFromDate((ARM_Date)fromDateStr);

    char toDateStr[11];
	Local_XLDATE2ARMDATE(toDate,toDateStr);
    toTime=mod.GetTimeFromDate((ARM_Date)toDateStr);

    startTime=mod.GetTimeFromDate(startDate);

    /*char endDateStr[11];
	Local_XLDATE2ARMDATE(endDate,endDateStr);*/
    endTime=mod.GetTimeFromDate(endDate);

    if( !(0.0<=fromTime && fromTime<=toTime && toTime<=startTime && startTime<=endTime) )
	{
		return ARM_KO;
	}

    return ARM_OK;
}

////////////////////////////////////////////
//// GP Vector Function
////////////////////////////////////////////
extern long ARMLOCAL_GP_Vector_Create(
	const VECTOR< double >& values,
	ARM_result&  result,
	long					objId )

{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& gpVector = NULL;
	try
	{			
		gpVector = new std::vector<double>(values);
		/// assign object
		if( !assignObject( gpVector, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	catch(Exception& x)
	{
		delete gpVector;
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}

////////////////////////////////////////////
//// GP Vector Function
////////////////////////////////////////////
extern long ARMLOCAL_GP_StrVector_Create(
	const VECTOR<string>& values,
	ARM_result&  result,
	long					objId )

{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GP_StrVector* gpVector = NULL;
	try
	{			
		gpVector = new ARM_GP_StrVector(values);
		/// assign object
		return( !assignObject( gpVector, result, objId )) ? ARM_KO : ARM_OK; 
	}
	catch(Exception& x)
	{
		delete gpVector;
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}


////////////////////////////////////////////
//// Trigonometric Matrix Function
////////////////////////////////////////////
extern long ARMLOCAL_TrigoMatrix(
	long n,
	double alpha,
	VECTOR< double >& matrix,
	long& nbRows,
	long& nbCols,
	ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GP_Matrix* gpMatrix = NULL;

	try
	{
		gpMatrix = TrigoMatrix(n, alpha);

		nbRows = gpMatrix->GetRowsNb();
		nbCols = gpMatrix->GetColsNb();

		matrix.resize(nbRows*nbCols);

		for (size_t i = 0; i < nbRows; ++i)
		{
			for (size_t j = 0; j < nbCols; ++j)
			{
				matrix[i*nbCols+j] = (*gpMatrix)(i,j);
			}
		}

		delete gpMatrix;
	}
	catch(Exception& x)
	{
		delete gpMatrix;
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;

}


////////////////////////////////////////////
//// Local covariance computation
////////////////////////////////////////////
extern long ARMLOCAL_LocalCovariance(
   long					modelId,
   string				underlyingType,
   double				fromDate,
   double				toDate,
   double				startDate1,
   string				endDate1Tenor,
   double				C_endDate1,
   long					endDate1Type,
   double				startDate2,
   string				endDate2Tenor,
   double				C_endDate2,
   long					endDate2Type,
   double				startDate3,
   string				endDate3Tenor,
   double				C_endDate3,
   long					endDate3Type,
   double				startDate4,
   string				endDate4Tenor,
   double				C_endDate4,
   long					endDate4Type,
   ARM_result&			result)
{
	/// Id checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_PricingModel* mod= NULL;
    try
	{
		ARM_PricingModelIR* oldModel = NULL;
		if( !GetObjectFromId( &mod, modelId, ARM_PRICINGMODEL) )
		{
			result.setMsg ("ARM_ERR: pricing model of a good type");
			return ARM_KO;
		}
       
		const char* calendar	= NULL;
		ARM_ZeroCurvePtr zeroCurve = mod->GetZeroCurve();
		
		if( zeroCurve != ARM_ZeroCurvePtr(NULL) )
			calendar = zeroCurve->GetCurrencyUnit()->GetCcyName();
		else
		{
			ARM_MultiAssetsModel* castedModel = dynamic_cast<ARM_MultiAssetsModel*> ( mod );
			if( castedModel )
				calendar = castedModel->GetRefModel()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
			else
			{
				result.setMsg ("ARM_ERR: Model has no ZeroCurve");
				return ARM_KO;
			}
		}
		
		char myStartDate1[20];
		char myStartDate2[20];
		char myStartDate3[20];
		char myStartDate4[20];
		
        // Convert dates
		/// the second argument can either be a date or a maturity
		Local_XLDATE2ARMDATE(startDate1, myStartDate1);
		Local_XLDATE2ARMDATE(startDate2, myStartDate2);
		Local_XLDATE2ARMDATE(startDate3, myStartDate3);
		Local_XLDATE2ARMDATE(startDate4, myStartDate4);
        ARM_Date startDate1(myStartDate1);
        ARM_Date startDate2(myStartDate2);
        ARM_Date startDate3(myStartDate3);
		ARM_Date startDate4(myStartDate4);

        ARM_Date endDate1(startDate1);
		ARM_GramFctorArg EndDate1Arg;
		if(endDate1Type == XL_TYPE_STRING)
		{
			EndDate1Arg = ARM_GramFctorArg(endDate1Tenor);
		}
		else
		{
			char myEndDate1[20];
			Local_XLDATE2ARMDATE(C_endDate1, myEndDate1);
            endDate1 = ARM_Date(myEndDate1);
			EndDate1Arg = ARM_GramFctorArg(endDate1);
		}

		
        ARM_Date endDate2(startDate2);
		ARM_GramFctorArg EndDate2Arg;
		if(endDate2Type == XL_TYPE_STRING)
		{
			EndDate2Arg = ARM_GramFctorArg(endDate2Tenor);
		}
		else
		{
			char myEndDate2[20];
			Local_XLDATE2ARMDATE(C_endDate2, myEndDate2);
            endDate2 = ARM_Date(myEndDate2);
			EndDate2Arg = ARM_GramFctorArg(endDate2);
		}

        ARM_Date endDate3(startDate3);
		ARM_GramFctorArg EndDate3Arg;
		if(endDate3Type == XL_TYPE_STRING)
		{
			EndDate3Arg = ARM_GramFctorArg(endDate3Tenor);
		}
		else
		{
			char myEndDate3[20];
			Local_XLDATE2ARMDATE(C_endDate3, myEndDate3);
            endDate3 = ARM_Date(myEndDate3);
			EndDate3Arg = ARM_GramFctorArg(endDate3);
		}

		ARM_Date endDate4(startDate4);
		ARM_GramFctorArg EndDate4Arg;
		if(endDate4Type == XL_TYPE_STRING)
		{
			EndDate4Arg = ARM_GramFctorArg(endDate4Tenor);
		}
		else
		{
			char myEndDate4[20];
			Local_XLDATE2ARMDATE(C_endDate4, myEndDate4);
            endDate4 = ARM_Date(myEndDate4);
			EndDate4Arg = ARM_GramFctorArg(endDate4);
		}


        if(endDate1Type == XL_TYPE_STRING || endDate1 != startDate1)
		    endDate1=ComputeEndDate( startDate1, EndDate1Arg, calendar );
        if(endDate2Type == XL_TYPE_STRING || endDate2 != startDate2)
		    endDate2=ComputeEndDate( startDate2, EndDate2Arg, calendar );
        if(endDate3Type == XL_TYPE_STRING || endDate3 != startDate3)
		    endDate3=ComputeEndDate( startDate3, EndDate3Arg, calendar );
		if(endDate4Type == XL_TYPE_STRING || endDate4 != startDate4)
		    endDate4=ComputeEndDate( startDate4, EndDate4Arg, calendar );

		double fromTime,toTime,startTime1,endTime1,startTime2,endTime2,startTime3,endTime3,startTime4,endTime4;
        if(ConvertDatesToTimes(*mod,fromDate,toDate,startDate1,endDate1,fromTime,toTime,startTime1,endTime1) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate1 <= endDate1 is expected");
			return ARM_KO;
		}
        if(ConvertDatesToTimes(*mod,fromDate,toDate,startDate2,endDate2,fromTime,toTime,startTime2,endTime2) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate2 <= endDate2 is expected");
			return ARM_KO;
		}
        if(ConvertDatesToTimes(*mod,fromDate,toDate,startDate3,endDate3,fromTime,toTime,startTime3,endTime3) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate3 <= endDate3 is expected");
			return ARM_KO;
		}
		if(ConvertDatesToTimes(*mod,fromDate,toDate,startDate4,endDate4,fromTime,toTime,startTime4,endTime4) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate4 <= endDate4 is expected");
			return ARM_KO;
		}

		double correl = mod->UnderlyingCovariance( underlyingType,
													fromTime,
													toTime,
													startTime1, 
													endTime1,
													startTime2,
													endTime2,
													startTime3,
													endTime3,
													startTime4,
													endTime4);

		result.setDouble(correl);
		return ARM_OK;

    }
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}	
}

////////////////////////////////////////////
//// Local covariance computation
////////////////////////////////////////////
extern long ARMLOCAL_LocalCorrelation(
   long					modelId,
   string				underlyingType,
   double				fromDate,
   double				toDate,
   double				startDate1,
   string				endDate1Tenor,
   double				C_endDate1,
   long					endDate1Type,
   double				startDate2,
   string				endDate2Tenor,
   double				C_endDate2,
   long					endDate2Type,
   double				startDate3,
   string				endDate3Tenor,
   double				C_endDate3,
   long					endDate3Type,
   double				startDate4,
   string				endDate4Tenor,
   double				C_endDate4,
   long					endDate4Type,
   ARM_result&			result)
{
	/// Id checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_PricingModel* mod= NULL;
    try
	{
		ARM_PricingModelIR* oldModel = NULL;
		if( !GetObjectFromId( &mod, modelId, ARM_PRICINGMODEL) )
		{
			result.setMsg ("ARM_ERR: pricing model of a good type");
			return ARM_KO;
		}
       
		const char* calendar	= NULL;
		ARM_ZeroCurvePtr zeroCurve = mod->GetZeroCurve();
		
		if( zeroCurve != ARM_ZeroCurvePtr(NULL) )
			calendar = zeroCurve->GetCurrencyUnit()->GetCcyName();
		else
		{
			ARM_MultiAssetsModel* castedModel = dynamic_cast<ARM_MultiAssetsModel*> ( mod );
			if( castedModel )
				calendar = castedModel->GetRefModel()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
			else
			{
				result.setMsg ("ARM_ERR: Model has no ZeroCurve");
				return ARM_KO;
			}
		}

        // Convert dates
		// the second argument can either be a date or a maturity
		char myStartDate1[20];
		char myStartDate2[20];
		char myStartDate3[20];
		char myStartDate4[20];
		
		Local_XLDATE2ARMDATE(startDate1, myStartDate1);
		Local_XLDATE2ARMDATE(startDate2, myStartDate2);
		Local_XLDATE2ARMDATE(startDate3, myStartDate3);
		Local_XLDATE2ARMDATE(startDate4, myStartDate4);

		ARM_GramFctorArg EndDate1Arg;
		if(endDate1Type == XL_TYPE_STRING)
		{
			EndDate1Arg = ARM_GramFctorArg(endDate1Tenor);
		}
		else
		{
			char myEndDate1[20];
			Local_XLDATE2ARMDATE(C_endDate1, myEndDate1);
			EndDate1Arg = ARM_GramFctorArg(ARM_Date(myEndDate1));
		}

		
		ARM_GramFctorArg EndDate2Arg;
		if(endDate2Type == XL_TYPE_STRING)
		{
			EndDate2Arg = ARM_GramFctorArg(endDate2Tenor);
		}
		else
		{
			char myEndDate2[20];
			Local_XLDATE2ARMDATE(C_endDate2, myEndDate2);
			EndDate2Arg = ARM_GramFctorArg(ARM_Date(myEndDate2));
		}

		ARM_GramFctorArg EndDate3Arg;
		if(endDate3Type == XL_TYPE_STRING)
		{
			EndDate3Arg = ARM_GramFctorArg(endDate3Tenor);
		}
		else
		{
			char myEndDate3[20];
			Local_XLDATE2ARMDATE(C_endDate3, myEndDate3);
			EndDate3Arg = ARM_GramFctorArg(ARM_Date(myEndDate3));
		}

		ARM_GramFctorArg EndDate4Arg;
		if(endDate4Type == XL_TYPE_STRING)
		{
			EndDate4Arg = ARM_GramFctorArg(endDate4Tenor);
		}
		else
		{
			char myEndDate4[20];
			Local_XLDATE2ARMDATE(C_endDate4, myEndDate4);
			EndDate4Arg = ARM_GramFctorArg(ARM_Date(myEndDate4));
		}

		ARM_Date endDate1= ComputeEndDate( ARM_Date(myStartDate1), EndDate1Arg, calendar );
		ARM_Date endDate2= ComputeEndDate( ARM_Date(myStartDate2), EndDate2Arg, calendar );
		ARM_Date endDate3= ComputeEndDate( ARM_Date(myStartDate3), EndDate3Arg, calendar );
		ARM_Date endDate4= ComputeEndDate( ARM_Date(myStartDate4), EndDate4Arg, calendar );

		double fromTime,toTime,startTime1,endTime1,startTime2,endTime2,startTime3,endTime3,startTime4,endTime4;
        if(ConvertDatesToTimes(*mod,fromDate,toDate,ARM_Date(myStartDate1),endDate1,fromTime,toTime,startTime1,endTime1) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate1 < endDate1 is expected");
			return ARM_KO;
		}
        if(ConvertDatesToTimes(*mod,fromDate,toDate,ARM_Date(myStartDate2),endDate2,fromTime,toTime,startTime2,endTime2) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate2 < endDate2 is expected");
			return ARM_KO;
		}
        if(ConvertDatesToTimes(*mod,fromDate,toDate,ARM_Date(myStartDate3),endDate3,fromTime,toTime,startTime3,endTime3) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate3 < endDate3 is expected");
			return ARM_KO;
		}
		if(ConvertDatesToTimes(*mod,fromDate,toDate,ARM_Date(myStartDate4),endDate4,fromTime,toTime,startTime4,endTime4) == ARM_KO)
		{
			result.setMsg ("ARM_ERR: fromDate <= toDate <= startDate4 < endDate4 is expected");
			return ARM_KO;
		}

		double correl = mod->UnderlyingCorrelation( underlyingType,
													fromTime,
													toTime,
													startTime1, 
													endTime1,
													startTime2,
													endTime2,
													startTime3,
													endTime3,
													startTime4,
													endTime4);

		result.setDouble(correl);
		return ARM_OK;

    }
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}




////////////////////////////////////////////
//// Function to create an event viewer representant
////////////////////////////////////////////
extern long ARMLOCAL_GetEventViewerRep(
	ARM_result& result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_EventViewerRep* eventViewer= NULL;

	try
	{
        eventViewer = new ARM_EventViewerRep;

		/// assign object
		if( !assignObject( eventViewer, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete eventViewer;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to reset message in an event viewer
////////////////////////////////////////////
extern long ARMLOCAL_GetEventViewer_ResetMssg(
	long eventViewerRepId,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_EventViewerRep* eventViewer = NULL;
		if( !GetObjectFromId( &eventViewer, eventViewerRepId, ARM_EVENT_VIEWER ) )
		{
			result.setMsg ("ARM_ERR: event viewer representant is not of a good type");
			return ARM_KO;
		};

		eventViewer->ResetMessage();

		result.setString( "Reseted Message" );
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to set the verbose mode in the event viewer
////////////////////////////////////////////
extern long ARMLOCAL_EventViewer_SetVerboseMode(
	bool verboseMode,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		string mssg( "Verbose Mode : " );
		mssg += ARM_EventViewerImp::SetVerboseMode( verboseMode ) ? "On": "Off";
		result.setString( mssg.c_str() );
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an error viewer representant
////////////////////////////////////////////
extern long ARMLOCAL_ErrViewer_Create(
	bool resetFlag,
	ARM_result& result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ErrViewer* errViewer  = NULL;

	try
	{
        errViewer = new ARM_ErrViewer(resetFlag);

		/// assign object
		if( !assignObject( errViewer, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete errViewer;
		x.DebugPrint();
		ARM_RESULT();
	}
}




//////////////////////////////////////////////////
//// Function to create a GP_Matrix
//////////////////////////////////////////////////

extern long ARMLOCAL_GP_Matrix_Create(
	const vector<double>& values,
	const long& rows,
	const long& cols,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GP_Matrix* matrix = NULL;
   
	try
	{ 
		matrix = new ARM_GP_Matrix( rows, cols );
		size_t i,j,k;
		
		for ( i=0, k=0;i<rows; ++i)
			for( j=0; j<cols; ++j )
				(*matrix)(i,j) = values[k++];
	
		/// assign object
		if( !assignObject( matrix, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete matrix ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to calculate a regression
//////////////////////////////////////////////////

extern long ARMLOCAL_GP_LeastSquareRegression(
	const VECTOR<double>& X,
	const long& rows,
	const long& cols,
	const VECTOR<double>& Y,
	VECTOR<double>& coeffs,
    ARM_result&				result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GP_Matrix XMatrix(rows,cols);
	std::vector<double> YVector(Y.size());
	std::vector<double>& CoeffsVector = NULL;
   
	try
	{ 
		size_t i,j,k;
		
		for ( i=0, k=0;i<rows; ++i)
			for( j=0; j<cols; ++j )
				XMatrix(i,j) = X[k++];

		for( i=0; i<Y.size(); ++i)
			YVector(i) = Y[i];

		CoeffsVector = LeastSquareRegression( XMatrix, YVector );

		coeffs.resize(CoeffsVector->size());

		for ( i=0;i<CoeffsVector->size();++i )
			coeffs[i] = (*CoeffsVector)(i);

		delete CoeffsVector;
	
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete CoeffsVector ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to calculate an acp
//////////////////////////////////////////////////

extern long ARMLOCAL_GP_ACP(
	const VECTOR<double>& InputMatrix,
	const long& rowsInput,
	const long& colsInput,
	VECTOR<double>& OuputMatrix,
	long& rowsOutput,
	long& colsOutput,
    ARM_result&				result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GP_Matrix matrix(rowsInput,colsInput);
	std::vector<double> eigenValues(colsInput);
	ARM_GP_Matrix* outMatrix = NULL;
   
	try
	{ 
		size_t i,j,k;
		
		for ( i=0, k=0;i<rowsInput; ++i)
			for( j=0; j<colsInput; ++j )
				matrix(i,j) = InputMatrix[k++];

		outMatrix = ACPTransformation( &matrix, eigenValues );

		rowsOutput = outMatrix->GetRowsNb()+1;
		colsOutput = outMatrix->GetColsNb();

		OuputMatrix.resize(rowsOutput*colsOutput);

		for ( i=0;i<outMatrix->GetColsNb();++i )
			OuputMatrix[i] = eigenValues(i);

		for ( i=0;i<outMatrix->GetRowsNb();++i )
			for ( j=0;j<outMatrix->GetColsNb();++j )
				OuputMatrix[(i+1)*outMatrix->GetColsNb()+j] = (*outMatrix)(i,j);

		delete outMatrix;
	
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete outMatrix ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_GP_Regression(
	const VECTOR<double>& Y,
	const VECTOR<double>& X,
	const long& nbRowsX,
	const long& nbColsX,
	const VECTOR<double>& XInter,
	const long& nbRowsXInter,
	const long& nbColsXInter,
	const string& RegMode,
	double span,
	VECTOR<double>& OuputVector,
    ARM_result& result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
   
	try
	{ 
		size_t i,j;
		
		ARM_GP_VectorPtr YVector(new std::vector<double>(Y.size()));
		for ( i=0;i<Y.size(); ++i)
			(*YVector)[i] = Y[i];

		ARM_GP_MatrixPtr XMatrix(new ARM_GP_Matrix(nbRowsX,nbColsX));
		for ( i=0;i<nbRowsX; ++i)
			for( j=0; j<nbColsX; ++j )
				(*XMatrix)(i,j) = X[i*nbColsX+j];

		ARM_GP_MatrixPtr XInterMatrix(new ARM_GP_Matrix(nbRowsXInter,nbColsXInter));
		for ( i=0;i<nbRowsXInter; ++i)
			for( j=0; j<nbColsXInter; ++j )
				(*XInterMatrix)(i,j) = XInter[i*nbColsXInter+j];

		ARM_Regression::RegressionMode regMode = (ARM_Regression::RegressionMode) ARM_ArgConv_RegMode.GetNumber(RegMode);
		ARM_RegressionPtr regression;
		if (regMode == ARM_Regression::LS)
		{
			regression = ARM_RegressionPtr(new ARM_LSRegression(YVector, XMatrix));
		}
		else if (regMode == ARM_Regression::LOESS)
		{
			regression = ARM_RegressionPtr(new ARM_LOESSRegression(YVector, XMatrix, span));
		}

		ARM_GP_VectorPtr res = regression->ComputeValues(XInterMatrix);

		OuputVector.resize(res->size());

		for ( i=0;i<OuputVector.size();++i )
			OuputVector[i] = (*res)[i];
	
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to calculate an integrated correlatoin
//////////////////////////////////////////////////

extern long ARMLOCAL_IntegratedCorrelation_Compute(
	const long& modelId,
	const double& tenor1,
	const vector<double>& tenors,
	const vector<double>& expiries,
	ARM_result& result,
	long objId  )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Surface* surface = NULL;

	try
	{
	    ARM_PricingModel* mod= NULL;
		if( !GetObjectFromId( &mod, modelId, ARM_PRICINGMODEL) )
		{
			result.setMsg ("ARM_ERR: pricing model of a good type");
			return ARM_KO;
		}

		ARM_GP_Matrix correlation( expiries.size(), tenors.size(), 0.0 );

		ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;
		for( size_t i=0; i<expiries.size(); ++i )
			for( size_t j=0; j<tenors.size(); ++j )
				correlation(i,j) = mod->UnderlyingCorrelation( "CMS", 
					0.0, 
					int( expiries[i]*365 ),
					int( expiries[i]*365 ), 
					int( (expiries[i]+tenor1)*365 ),  
					int( expiries[i]*365 ), 
					int( (expiries[i]+tenors[i])*365 ),
					int( expiries[i]*365 ), 
					int( (expiries[i]+tenor1)*365 ),
					int( expiries[i]*365 ), 
					int( (expiries[i]+tenor1)*365 ) );  
		surface = new ARM_SurfaceWithInterpol( std::vector<double>(expiries), std::vector<double>(tenors), correlation,type);

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


//////////////////////////////////////////////////
//// Function to calculate an integrated correlatoin
//////////////////////////////////////////////////

extern long ARMLOCAL_FxMktToTotemCalibrate(
	const long& portfolioId,
	const double& atmVol,
	const vector<double>& deltaCalls,
	const vector<double>& deltasPuts,
	const long& modelId,
	const long& initpointId,
	const long& lowBoundId,
	const long& upBoundId,
	const double& maxIter,
	const string& algoTypeStr,
	const double& tolerance,
	const string& OutPutId,
	ARM_result& result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* Res = NULL;

	ARM_FXMktDataToTotemFormat* fxMkt = NULL;
	ARM_AutoCleaner< ARM_FXMktDataToTotemFormat > HoldFXMKT(fxMkt);

	ARM_BSModel* bsModel = NULL;
    ARM_AutoCleaner< ARM_BSModel > HoldBSM(bsModel);

	ARM_StdPortfolio*  portfolio = NULL;
	ARM_AutoCleaner< ARM_StdPortfolio > HoldPF(portfolio);

	std::vector<double>&  initVariable = NULL;
	ARM_AutoCleaner< std::vector<double> > HoldIV(initVariable);

	std::vector<double>&  lowBound = NULL;
	ARM_AutoCleaner< std::vector<double> > HoldLB(lowBound);

	std::vector<double>&  upBound = NULL;
	ARM_AutoCleaner< std::vector<double> > HoldUB(upBound);

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Fx Mkt Data to Totem Data" );

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &bsModel, modelId ,"black & Scholes", result ) ) return ARM_KO;

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &portfolio, portfolioId ,"Portfolios", result ) ) return ARM_KO;

		//Initial Point
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &initVariable,initpointId,"init variable",result ) ) return ARM_KO;

		//fundLevrage
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &lowBound, lowBoundId,"lower Bound ",result ) ) return ARM_KO;

		//fundLevrage
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &upBound,upBoundId,"upper Bound",result ) ) return ARM_KO;

		int algoType = ARM_ArgConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag.GetNumber(algoTypeStr);

		std::vector<double> vdeltaCalls(deltaCalls);
		std::vector<double> vdeltaPuts(deltasPuts);

		ARM_FXMktDataToTotemFormat fxMkt(*portfolio,
			vdeltaCalls,
			vdeltaPuts,
			atmVol,
			*bsModel,
			*initVariable,
			*lowBound,
			*upBound,
			maxIter,
			algoType,
			tolerance);

		//Res = (std::vector<double>&)fxMkt.GetUnKnown().Clone();
		if(OutPutId == string("MKTMODEL"))
			Res = (ARM_BSModel*)fxMkt.GetMktModel()->Clone();
		
		else if(OutPutId == string("GENCURVE")){
			std::vector<double> values(fxMkt.GetUnKnown());
			values/=100.0;
			std::vector<double> deltas(values.size());
			for(size_t i(0); i<values.size(); ++i)
				deltas[i] = i;
			Res = new ARM_Curve(deltas,values, new ARM_LinInterpCstExtrapolDble );
		}
		else{
			result.setMsg ("ARM_ERR: OutPut id: string invalid");
			return ARM_KO;
		}
		/// assign object
		if( !assignObject( Res, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete Res;
		x.DebugPrint();
		ARM_RESULT();
	}
}