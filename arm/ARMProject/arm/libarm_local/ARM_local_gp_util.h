/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_local_gp_util.h
 *
 *  \brief file for the various utils in the generic pricer
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date January 2004
 */

#ifndef ARMLOCAL_GP_UTIL_H
#define ARMLOCAL_GP_UTIL_H

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>

class ARM_GramFctorArg;

class ARM_PricingModel;

////////////////////////////////////////////
//// XL Dates to GP Times convertion
////////////////////////////////////////////
long ConvertDatesToTimes(
	ARM_PricingModel& mod,
	double       fromDate,
	double       toDate,
	double       startDate,
	double       endDate,
	double&      fromTime,
	double&      toTime,
	double&      startTime,
	double&      endTime);

////////////////////////////////////////////
//// function to create a gp vector
////////////////////////////////////////////
extern long ARMLOCAL_GP_Vector_Create(
	const vector<double>& values,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_GP_StrVector_Create(
	const VECTOR<string>& values,
	ARM_result&  result,
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Trigonometric Function
////////////////////////////////////////////
extern long ARMLOCAL_TrigoMatrix(
	long n,
	double alpha,
	VECTOR< double >& matrix,
	long& nbRows,
	long& nbCols,
	ARM_result&  result);


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
   ARM_result&			result);


////////////////////////////////////////////
//// Local correlation computation
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
   ARM_result&			result);


////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
extern long ARMLOCAL_GetEventViewerRep(
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// function to reset message in the event viewer
////////////////////////////////////////////
extern long ARMLOCAL_GetEventViewer_ResetMssg(
	long eventViewerRepId,
	ARM_result& result );

////////////////////////////////////////////
//// function to set verbose mode in the event viewer
////////////////////////////////////////////
extern long ARMLOCAL_EventViewer_SetVerboseMode(
	bool verboseMode,
	ARM_result& result );

////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
extern long ARMLOCAL_ErrViewer_Create(
	bool resetFlag,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID);



////////////////////////////////////////////
//// function to create a gp matrix
////////////////////////////////////////////
extern long ARMLOCAL_GP_Matrix_Create(
	const vector<double>& values,
	const long& rows,
	const long& cols,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// function to compute a regression
////////////////////////////////////////////
extern long ARMLOCAL_GP_LeastSquareRegression(
	const VECTOR<double>& X,
	const long& rows,
	const long& cols,
	const VECTOR<double>& Y,
	VECTOR<double>& coeffs,
    ARM_result&				result);

////////////////////////////////////////////
//// function to compute an ACP
////////////////////////////////////////////
extern long ARMLOCAL_GP_ACP(
	const VECTOR<double>& InputMatrix,
	const long& rowsInput,
	const long& colsInput,
	VECTOR<double>& OuputMatrix,
	long& rowsOutput,
	long& colsOutput,
    ARM_result&				result);


////////////////////////////////////////////
//// function to compute a Regression with
//// LS or LOESS algorithm
////////////////////////////////////////////
extern long ARMLOCAL_GP_Regression(
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
    ARM_result& result);


////////////////////////////////////////////
//// function to compute the integrated correlation
////////////////////////////////////////////
extern long ARMLOCAL_IntegratedCorrelation_Compute(
	const long& modelId,
	const double& tenor1,
	const vector<double>& tenors,
	const vector<double>& expiries,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// function to convert Fx Mkt data to TOTEM format
////////////////////////////////////////////
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
	long objId = ARM_NULL_OBJECT_ID );

#endif
