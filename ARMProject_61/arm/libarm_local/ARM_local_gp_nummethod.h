/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_nummethod.h,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_local_gp_nummethod.h
 *
 *  \brief file for the numerical method part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef ARMLOCAL_GP_NUMMETHOD_H
#define ARMLOCAL_GP_NUMMETHOD_H

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_result.h"

extern string ModelNumMethodToClass(long modelId);

////////////////////////////////////////////
//// Function to create a backward induction
///		numerical method
////////////////////////////////////////////
extern long ARMLOCAL_BINumMethod_Create(
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


//// Function to create a forward induction
///		numerical method
////////////////////////////////////////////
extern long ARMLOCAL_FINumMethod_Create(
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);

//// Function to create a mixte induction
///		numerical method
////////////////////////////////////////////
extern long ARMLOCAL_MixteNumMethod_Create(
	ARM_result&	result, 
	long objId );


////////////////////////////////////////////
//// Function to set a numerical method to a model
///		numerical method
////////////////////////////////////////////
extern long ARMLOCAL_SetNumMethod(
	long modelId,
	long numMethodId,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to set numeraire to model
////////////////////////////////////////////
extern long ARMLOCAL_SetNumeraire(
	long modelId,
	long numMethodId,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a Numeraire
////////////////////////////////////////////
extern long ARMLOCAL_Numeraire_Create(
    long                    numeraireType,
    const VECTOR<double>&   numeraireTimes,
	ARM_result&	            result, 
	long                    objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a Tree method
////////////////////////////////////////////
extern long ARMLOCAL_TreeMethod_Create(
    int         nbSteps,
    double      stdDev,
    double      minStdDev,
    int         nbMinSteps,
    bool        isTree1GForced,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a Monte Carlo
////////////////////////////////////////////
extern long ARMLOCAL_MCMethod_Create(
    const size_t&			itersNb,
    const int&				fixStep,
	const VECTOR<long>&		randGensId,
	const long&				C_SamplerType,
	const VECTOR<double>&	C_SamplerDatas,
	const long&				C_SchedulerType,
    const VECTOR<double>&	C_SchedulerDatas,
	const long&				C_ImpSamplerType,
	const VECTOR<double>&	C_ImpSamplerDatas,
	const long&				C_PathSchemeType,
	const long&				ExercBoundCalcId,
	const size_t&			MaxBucketSize,
	ARM_result&				result,
	long			objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a PDE
////////////////////////////////////////////
//Initial version
extern long ARMLOCAL_PDEMethod_Create(
	const string&	MethodName,				// Method Name
    const vector<double>& SchedulerData,	// Scheduler Data
	const int&	SchedulerDataNbRows,			// SchedulerNbRows
	const int&	SchedulerDataNbCols,			// SchedulerNbCols
	const int&	SpaceItersNb,			// Space Iter Nb 
	const string&	GridType,				// Grid Type
	const vector<double>& GridData,			// Grid Data
	const int&	GridDataNbRows,			// GridNbRows
	const int&	GridDataNbCols,			// GridNbCols
	const int&	YGridItersNb,			// YGridIterNb
	const int&	ZGridItersNb,			// ZGridIterNb
	const double&	Theta1,					// Theta1
	const double&	Theta2,					// Theta2
	const double&	Theta3,					// Theta3
	const string&	BoundaryConditionName,	// Bound Cond
	const double&	lambda,					// Lambda
	ARM_result&		result,
	long			objId = ARM_NULL_OBJECT_ID);

//ND generalisation
extern long ARMLOCAL_PdeND_Create(
	string			C_MethodName,
    const long&     C_SchedulerType,
	long			C_SchulerData,
	long			C_SpaceData,
	long			C_SchemeData,
	string			C_BoundCondName,
	ARM_result&		result, 
	long			objId  = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a random nb generator
////////////////////////////////////////////
extern long ARMLOCAL_RandGen_Create(
	const string& genMode,
	const string& algo,
	const long& baseGen1Id,
	const long& baseGen2Id,
	const int& seed,
	const int& dim,
	const int& factorDim,
	const int& nbOfPoints,
	const double& nbStdDevs,
	const int& firstNbTimes,
	const int& firstNbDims,
	const string& order,
	const int& firstSimulations,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a simple random nb generator
////////////////////////////////////////////
extern long ARMLOCAL_SimpleRandGen_Create(
	const string& genType1,
	const string& genType2,
	const string& algo1,
	const string& algo2,
	const int& firstNbTimes,
	const int& firstNbDims,
	const string& isAntithetic,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to draw nbs from a random nb generator
////////////////////////////////////////////
extern long ARMLOCAL_RandGen_DrawVector(
	long randomGenId,
	int size,
	VECTOR<double>& data,
	ARM_result&	result );


///////////////////////////////////////////////
//// Get the numerical method of a GP model
///////////////////////////////////////////////
extern long ARMLOCAL_GetNumMethodFromModel(
        long modelId,
        ARM_result&	result, 
        long        objId= ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Set the spot proba computation flag of a tree method
////////////////////////////////////////////
extern long ARMLOCAL_Local_Tree_SetProbaFlag(
	long treeId,
	bool isSpotProba,
	ARM_result&	result );


////////////////////////////////////////////
//// Function to create an Andersen
////////////////////////////////////////////
extern long ARMLOCAL_AMCAndersen_Create(
	const double& ItersNb, 
	const bool& sortedMaximization,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a Longstaff&Schwartz
////////////////////////////////////////////
extern long ARMLOCAL_AMCLongstaffSchwartz_Create(
	const double& ItersNb,
	const string& RegMode,
	const double& Span,
	const string& isAutomatic,
	const int& degree,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a Tree 1D method
////////////////////////////////////////////
extern long ARMLOCAL_TreeND_Create(
         const long&                C_NbDims,
         const long&                C_SchedulerType,
         const VECTOR<double>&      C_SchedulerDatas,
         const long&                C_SamplerType,
         const VECTOR<double>&      C_SamplerDatas,
         const long&                C_TruncatorType,
         const VECTOR<double>&      C_TruncatorDatas,
         const bool&                C_ProbasFlag,
         const long&                C_ReconnectorType,
         const long&                C_SmootherType,
         ARM_result&	result, 
         long        objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_CFMethod_Create(
		 const string&	MethodeName,
		 const long&	matrixId,	
         ARM_result&	result, 
         long        objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_ImpSampler_Optimize(
         const long&                C_GenSecId,
         const long&                C_ModelId,
         const long&                C_InitGuessId,
		 const double&				C_InitGuess,
         const long&                C_LowerBoundId,
		 const double&              C_LowerBound,
		 const long&                C_UpperBoundId,
		 const double&              C_UpperBound,
		 const string&				C_withMC,
		 const long&				C_nbSteps,
		 const string&				C_bootstrap,
         ARM_result&	result, 
         long        objId = ARM_NULL_OBJECT_ID );


#endif