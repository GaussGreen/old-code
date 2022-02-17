/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file typedef.h
 *  \brief global typedef of namespace ARM
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPNUMLIB_TYPEDEF_H
#define _INGPNUMLIB_TYPEDEF_H

#include "gpbase/port.h"
#include "gpbase/countedptr.h"
#include "gpbase/functor.h"
#include "gpbase/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// struct forward declaration
struct ARM_ManipulatorMethod;

class ARM_RandomGenerator;
class ARM_MomentFunc;
class ARM_Regression;

/// reference counted pointor
typedef	ARM_CountedPtr< ARM_RandomGenerator >			ARM_RandomGeneratorPtr;
typedef	ARM_CountedPtr< ARM_ManipulatorMethod >			ARM_ManipulatorMethodPtr;
typedef ARM_CountedPtr< ARM_MomentFunc>					ARM_MomentFuncPtr;
typedef	vector< ARM_RandomGeneratorPtr >				ARM_RandomGeneratorPtrVector;
typedef vector < ARM_MomentFuncPtr >					ARM_MomentFuncPtrVector;


/// forward declaration in ARM namespace
template <typename Func, typename T=double>       class T_SolverWithInitialGuess;
class FunctionToSolveWithDerivative;
///	typedef BootStrap1DSolverType SolverType;
typedef T_SolverWithInitialGuess<FunctionToSolveWithDerivative> ModifiedNRSolver;
typedef ARM_CountedPtr< ModifiedNRSolver > ModifiedNRSolverPtr;
typedef ARM_CountedPtr< ARM_Regression > ARM_RegressionPtr;

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
