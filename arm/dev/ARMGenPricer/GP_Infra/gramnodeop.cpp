/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramnodeop.cpp,v $
 * Revision 1.1  2003/11/23 16:41:37  ebenhamou
 * Initial revision
 *
 *
 */


/*! \file gramnodeop.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/gramnodeop.h"


CC_BEGIN_NAMESPACE( ARM )

/// negate
CC_NS( ARM_GP, negate<double> ) UnaryOpNegate;

/// logical_not
CC_NS( ARM_GP, logical_not<double> ) UnaryOpLogicalNot;

/// plus
CC_NS( ARM_GP, plus<double> ) BinaryOpPlus;

/// minus
CC_NS( ARM_GP, minus<double> ) BinaryOpMinus;

/// multiplies
CC_NS( ARM_GP, multiplies<double> ) BinaryOpMultiplies;

/// divides
CC_NS( ARM_GP, divides<double> ) BinaryOpDivides;

/// modulus
CC_NS( ARM_GP, modulus<double> ) BinaryOpModulus;

/// equal_to
CC_NS( ARM_GP, equal_to<double> ) BinaryOpEqualTo;

/// not_equal_to
CC_NS( ARM_GP, not_equal_to<double> ) BinaryOpNotEqualTo;

/// less
CC_NS( ARM_GP, less<double> ) BinaryOpLess;

/// greater
CC_NS( ARM_GP, greater<double> ) BinaryOpGreater;

/// less_equal
CC_NS( ARM_GP, less_equal<double> ) BinaryOpLessEqual;

/// greater_equal
CC_NS( ARM_GP, greater_equal<double> ) BinaryOpGreaterEqual;

/// logical_and
CC_NS( ARM_GP, logical_and<double> ) BinaryOpLogicalAnd;

/// logical_or
CC_NS( ARM_GP, logical_or<double> ) BinaryOpLogicalOr;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


