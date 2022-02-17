/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctor.h,v $
 * Revision 1.1  2003/15/09 18:53:24  ebenhamou
 *  filtering of blank
 *
 *
 */


/*! \file gramfunctorconv.h
 *
 *  \brief this files deals with all the convention for the 
 *		expression node class of the grammar
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPINFRA_GRAMFUNCTORCONV_H
#define _INGPINFRA_GRAMFUNCTORCONV_H

#include "gpbase/port.h"
#include "gramfunctorarg.h"
#include "dates.h"



/// forward declaration 
class ARM_Currency;


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_ExpNode;
class ARM_PricingModel;


/// gets and checks evaluation time
double GetAndCheckEvalTime( double evalDate, ARM_PricingModel* mod, const string& funcName );

/// function to compute a reset or pay date
ARM_Date ComputeResetOrPayDate(
	const ARM_Date& advDate,
	const ARM_Date& arrDate,
    int timing,
	char* calendar,
    const ARM_GramFctorArg& argGapOrDate,
    ARM_PricingModel* mod, 
	const string& curveName,
    const string& funcName);

/// function to compute a settlement date
ARM_Date ComputeSettlementDate( 
	const ARM_Date& startDate, 
	char* calendar, 
	const ARM_GramFctorArg& argGapOrDate,
	int defaultGap,
    const string& funcName );

/// function to set the settlement date
void SetSettlementDate( 
	vector< ARM_ExpNodePtr >& nodes, /// changed potentially
	const ARM_Date& startDate, 
	const string& calendar,
	int defaultGap,
	ARM_GramFctorArgVector& args, /// changed potentially
	size_t i,
    const string& funcName );

/// function to compute the index term
string ComputeIndexTerm( const ARM_GramFctorArg& argTerm, const ARM_GramFctorArg& argStart,
                         const ARM_GramFctorArg& argEnd,const string& funcName);

/// function to set the index term
void SetIndexTerm( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i,
                   size_t startIdx, size_t endIdx, ARM_PricingModel* mod, const string& curveName,
                   const string& funcName);

/// compute the reset days gap getting either the default or the one given
double GetResetDaysGap( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& CurveName );

/// Set reset days gap
void SetResetDaysGap( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName );

/// sets the end date
void SetEndDateResetCal( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName );
void SetEndDatePayCal( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName );

/// computes the end date
ARM_Date ComputeEndDate( const ARM_Date& startDate, const ARM_GramFctorArg& arg,const char* calendar );

/// Very simple structures to handle fixed/floating DC method
struct GetDayCountFtor
{
    GetDayCountFtor::GetDayCountFtor(ARM_Currency* ccy) : itsCcy(ccy) {}
    virtual long operator() () const = 0;
    ARM_Currency* itsCcy;
};

struct GetFixedDayCountFtor : public GetDayCountFtor
{
    GetFixedDayCountFtor::GetFixedDayCountFtor(ARM_Currency* ccy)
        : GetDayCountFtor(ccy) {}
    virtual long operator() () const;
};

struct GetFloatDayCountFtor : public GetDayCountFtor
{
    GetFloatDayCountFtor::GetFloatDayCountFtor(ARM_Currency* ccy)
        : GetDayCountFtor(ccy) {}
    virtual long operator() () const;
};

/// get the day count of a leg
int GetDayCount( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& CurveName, const GetDayCountFtor& getDayCount);

/// set method for day count of a leg
bool SetDayCount(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, const GetDayCountFtor& getDayCount);

/// Very simple structures to handle fixed/floating frequency method
struct GetFrequencyFtor
{
    GetFrequencyFtor::GetFrequencyFtor(ARM_Currency* ccy) : itsCcy(ccy) {}
    virtual double operator() () const = 0;
    ARM_Currency* itsCcy;
};

struct GetFixedFrequencyFtor : public GetFrequencyFtor
{
    GetFixedFrequencyFtor::GetFixedFrequencyFtor(ARM_Currency* ccy)
        : GetFrequencyFtor(ccy) {}
    virtual double operator() () const;
};

struct GetFloatFrequencyFtor : public GetFrequencyFtor
{
    GetFloatFrequencyFtor::GetFloatFrequencyFtor(ARM_Currency* ccy)
        : GetFrequencyFtor(ccy) {}
    virtual double operator() () const;
};


// get term of a leg
double GetTerm(const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& CurveName, const string funcName, const GetFrequencyFtor& getFrequency);

// set term of a leg
bool SetTerm(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, const GetFrequencyFtor& getFrequency);

// set term of a leg
bool SetTerm(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, double term);

/// get frequency of a leg
int GetFrequency( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& CurveName, const string funcName, const GetFrequencyFtor& getFrequency);

/// set method for frequency of a leg with a functor
bool SetFrequency(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, const GetFrequencyFtor& getFrequency);

/// set method for frequency of a leg with a frequency value
bool SetFrequency(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, int frequency);

/// set method for frequency of a leg with a frequency value (string)
bool SetFrequencyToGramFctorArg(vector< ARM_ExpNodePtr >& nodes, ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, const string& freqStr);

ARM_GP_CurvePtr GetCurve(const ARM_GramFctorArg& arg);

/// Convert a double, vector or curve to a vector w.r.t. times vector (and an offset in case of input vector)
ARM_VectorPtr GetProfileValues(const ARM_GramFctorArg& arg, size_t offset,
                               const ARM_GP_Vector& times, const string& msg);

CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

