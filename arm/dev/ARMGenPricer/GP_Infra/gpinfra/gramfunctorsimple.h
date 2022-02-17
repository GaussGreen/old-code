/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctorsimple.h,v $
 * Revision 1.1  2003/15/09 18:53:24  ebenhamou
 *  filtering of blank
 *
 *
 */


/*! \file gramfunctorsimple.h
 *
 *  \brief contains all very simple gramfunctor
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_GRAMFUNCTORSIMPLE_H
#define _INGPINFRA_GRAMFUNCTORSIMPLE_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gramfunctorbase.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// general class for unary functor
struct ARM_GP_UnaryOp : public ARM_GramFctor
{
	ARM_GP_UnaryOp( const pDbleUnaryFunc& Op, const string& FuncName ): ARM_GramFctor(), itsOp(Op), itsFuncName( FuncName ) {}
    ARM_GP_UnaryOp(const ARM_GP_UnaryOp& rhs) : ARM_GramFctor(), itsOp(rhs.itsOp), itsFuncName( rhs.itsFuncName ) {}
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_UnaryOp( *this ) ); }
private:
	ARM_GramFctorArg OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	ARM_GramFctorArg OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

	pDbleUnaryFunc itsOp;
	/// name for error writting
	string itsFuncName;
};

/// general class for binary functor
struct ARM_GP_BinaryOpWithDates : public ARM_GramFctor
{
	ARM_GP_BinaryOpWithDates( const pDbleBinaryFunc& Op, const string& FuncName ): ARM_GramFctor(), itsOp(Op), itsFuncName( FuncName ) {}
    ARM_GP_BinaryOpWithDates(const ARM_GP_BinaryOpWithDates& rhs) : ARM_GramFctor(), itsOp(rhs.itsOp), itsFuncName( rhs.itsFuncName ) {}
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_BinaryOpWithDates( *this ) ); }
	virtual const string& GetFuncName() { return itsFuncName; }
private:
	ARM_GramFctorArg OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	ARM_GramFctorArg OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

	pDbleBinaryFunc itsOp;
	/// name for error writting
	string itsFuncName;
};


/// Power operator: a little specific hence the complete definition
struct ARM_GP_PowVector : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_PowVector() : ARM_GramFctor() {}
    ARM_GP_PowVector(const ARM_GP_PowVector& rhs) : ARM_GramFctor(rhs) {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_PowVector( *this ) ); }
private:
	ARM_GramFctorArg OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	ARM_GramFctorArg OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    static string itsFuncName;
};

/// If operator: a little specific hence the complete definition
struct ARM_GP_IfVector : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_IfVector() : ARM_GramFctor() {}
    ARM_GP_IfVector(const ARM_GP_IfVector& rhs) : ARM_GramFctor(rhs) {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_IfVector( *this ) ); }
private:
	ARM_GramFctorArg OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	ARM_GramFctorArg OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    static string itsFuncName;
};

/// PV operator
struct ARM_GP_PV : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_PV() : ARM_GramFctor() {}
    ARM_GP_PV(const ARM_GP_PV& rhs) : ARM_GramFctor(rhs) {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_PV( *this ) ); }
private:
    static string itsFuncName;
};

/// UnPay operator
struct ARM_GP_UnPay : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_UnPay() : ARM_GramFctor() {}
    ARM_GP_UnPay(const ARM_GP_UnPay& rhs) : ARM_GramFctor(rhs) {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_UnPay( *this ) ); }
private:
	/// name for error writting
    static string itsFuncName;
};


/// Day Count Fraction
struct ARM_GP_DCF : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_DCF() : ARM_GramFctor() {}
    ARM_GP_DCF(const ARM_GP_DCF& rhs) : ARM_GramFctor(rhs) {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_DCF( *this ) ); }
private:
    static string itsFuncName;
};


/// Exercise operator
struct ARM_GP_Exercise : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_Exercise();
    ARM_GP_Exercise(const ARM_GP_Exercise& rhs);
	~ARM_GP_Exercise();

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Exercise( *this ) ); }
	inline ARM_ExerciseBoundary * GetExerciseBoundary( void ) { return itsExerciseBoundary; }
private:
    static string itsFuncName;
	ARM_ExerciseBoundary * itsExerciseBoundary;
	ARM_GP_MatrixPtr itsMat;
};

/// Frontier operator
struct ARM_GP_Frontier : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_Frontier();
    ARM_GP_Frontier(const ARM_GP_Frontier& rhs);
	~ARM_GP_Frontier();

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Frontier( *this ) ); }
private:
    static string itsFuncName;
};


/// Trigger operator
struct ARM_GP_Trigger : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_Trigger();
    ARM_GP_Trigger(const ARM_GP_Trigger& rhs);
	~ARM_GP_Trigger();

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Trigger( *this ) ); }
private:
    static string itsFuncName;
	/// contains testVar2, so that it is computed in fwd loop only
	ARM_GP_VectorPtr testVar2;
};

struct ARM_GP_MFactorFctor : public ARM_GramFctor
{
	ARM_GP_MFactorFctor();
	ARM_GP_MFactorFctor( const ARM_GP_MFactorFctor& rhs );
	~ARM_GP_MFactorFctor();

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_MFactorFctor( *this ) ); }

private:
    static string itsFuncName;
};

struct ARM_GP_SumSerieVector : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_SumSerieVector() : ARM_GramFctor() {}
    ARM_GP_SumSerieVector(const ARM_GP_SumSerieVector& rhs) : ARM_GramFctor(rhs) {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_SumSerieVector( *this ) ); }
private:
	ARM_GramFctorArg OperatorInPlace( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
	ARM_GramFctorArg OperatorWithClone( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

    static string itsFuncName;
};

CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

