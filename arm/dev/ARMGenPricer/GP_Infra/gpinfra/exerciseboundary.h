/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file exerciseboundary.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPINFRA_EXERCISEBOUNDARY_H
#define _INGPINFRA_EXERCISEBOUNDARY_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/assignop.h"

#include "gpnumlib/typedef.h"
#include "gpnumlib/regression.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ExerciseBoundary : public ARM_RootObject
{
private:
	int itsIsAutomatic;
	int itsDegree;
public:

	ARM_ExerciseBoundary(int isAutomatic, int degree)
		:
		itsIsAutomatic(isAutomatic),
		itsDegree(degree) 
	{};

	//virtual ~ARM_ExerciseBoundary {};

	ARM_ExerciseBoundary(const ARM_ExerciseBoundary&  rhs )
		:
		itsIsAutomatic(rhs.itsIsAutomatic),
		itsDegree(rhs.itsDegree) 
	{};

	/// Node manipulation
	virtual void EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector ) = NULL;

	/// Accessor to Exercise Boundary
	virtual ARM_GP_VectorPtr GetExerciseBoundary( void ) const = 0;

	virtual int IsAutomatic() const  { return itsIsAutomatic; };
	virtual int Degree() const { return itsDegree; };

};

class ARM_AndersenExerciseBoundary : public ARM_ExerciseBoundary
{
private: 
	double itsValue;
	bool itsExerciseIfPayoffIsGreaterThanValue;

public: 
	/// Node manipulation
	virtual void EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// Constructors/Destructor
	ARM_AndersenExerciseBoundary( double value = 0.0, bool exerciseType = true ) 
		:
		ARM_ExerciseBoundary(0,0),
		itsValue(value),
		itsExerciseIfPayoffIsGreaterThanValue(exerciseType)
	{};

	virtual ~ARM_AndersenExerciseBoundary() {};

	ASSIGN_OPERATOR(ARM_AndersenExerciseBoundary)

	ARM_AndersenExerciseBoundary(const ARM_AndersenExerciseBoundary& rhs)
		:
		ARM_ExerciseBoundary(rhs),
		itsValue(rhs.itsValue),
		itsExerciseIfPayoffIsGreaterThanValue(rhs.itsExerciseIfPayoffIsGreaterThanValue)
	{};

	/// Accessors
	inline double getValue() const { return itsValue; };
	void setValue( const double & value ) { itsValue = value; };
	virtual ARM_GP_VectorPtr GetExerciseBoundary( void ) const;
	inline bool ExerciseIfPayoffIsGreaterThanValue() const { return itsExerciseIfPayoffIsGreaterThanValue; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_AndersenExerciseBoundary(*this); };
    virtual string toString(const string& indent="", const string& nextIndent="") const; 
};


class ARM_LSExerciseBoundary : public ARM_ExerciseBoundary
{
private: 
	void CopyNoCleanUp(const ARM_LSExerciseBoundary& rhs);
    void CleanUp();
	ARM_RegressionPtr itsRegression;

	ARM_GP_VectorPtr ComputePseudoContinuationValues( const ARM_GP_MatrixPtr& StatesVector );

public: 
	/// Node manipulation
	virtual void EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// Constructors/Destructor
	ARM_LSExerciseBoundary( const ARM_RegressionPtr& Regression, int isAutomatic, int degree ) 
		:
		ARM_ExerciseBoundary(isAutomatic, degree),
		itsRegression(Regression)
	{};

	~ARM_LSExerciseBoundary() {};

	ASSIGN_OPERATOR(ARM_LSExerciseBoundary)

	ARM_LSExerciseBoundary( const ARM_LSExerciseBoundary& rhs ) 
		:
		ARM_ExerciseBoundary(rhs),
		itsRegression(rhs.itsRegression)
	{};


	/// Accessors
	virtual ARM_GP_VectorPtr GetExerciseBoundary( void ) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_LSExerciseBoundary(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const; 
};


class ARM_SmoothTreeExerciseBoundary : public ARM_ExerciseBoundary
{
private:
    /// Values for payoff smoothing
	ARM_GP_VectorPtr itsSmoothValues;

    /// Interpolated exercise boundary (not computed at the moment)
    ARM_GP_VectorPtr itsExerValues;

public: 
	/// Node manipulation
	virtual void EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// Constructors/Destructor
    ARM_SmoothTreeExerciseBoundary( const ARM_GP_VectorPtr& smoothValues, const ARM_GP_VectorPtr& exerValues) 
		: 
		ARM_ExerciseBoundary(0,0),
		itsSmoothValues(smoothValues), 
		itsExerValues(exerValues)
	{};

	virtual ~ARM_SmoothTreeExerciseBoundary()
	{
	}

	ASSIGN_OPERATOR(ARM_SmoothTreeExerciseBoundary)

	ARM_SmoothTreeExerciseBoundary( const ARM_SmoothTreeExerciseBoundary& rhs) 
		: 
		ARM_ExerciseBoundary(rhs),
		itsSmoothValues(rhs.itsSmoothValues), 
		itsExerValues(rhs.itsExerValues)
	{};

	/// Accessors
    virtual ARM_GP_VectorPtr GetExerciseBoundary( void ) const { return itsExerValues; }

	/// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_SmoothTreeExerciseBoundary(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_SmoothTreeExerciseBoundary"; }
};

CC_END_NAMESPACE()

#endif
