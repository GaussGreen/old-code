/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramnodeop.h,v $
 * Revision 1.1  2003/10/22 13:15:24  ebenhamou
 * Initial version
 *
 *
 */


/*! \file gramnodeop.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_GRAMNODEOP_H
#define _INGPINFRA_GRAMNODEOP_H

#include "gpbase/port.h"
#include "functordef.h"

CC_BEGIN_NAMESPACE( ARM )

/// negate
extern CC_NS( ARM_GP, negate )<double> UnaryOpNegate;

/// logical_not
extern CC_NS( ARM_GP, logical_not )<double> UnaryOpLogicalNot;

/// plus
extern CC_NS( ARM_GP, plus )<double> BinaryOpPlus;

/// minus
extern CC_NS( ARM_GP, minus )<double> BinaryOpMinus;

/// multiplies
extern CC_NS( ARM_GP, multiplies )<double> BinaryOpMultiplies;

/// divides
extern CC_NS( ARM_GP, divides )<double> BinaryOpDivides;

/// modulus
extern CC_NS( ARM_GP, modulus )<double> BinaryOpModulus;

/// equal_to
extern CC_NS( ARM_GP, equal_to )<double> BinaryOpEqualTo;

/// not_equal_to
extern CC_NS( ARM_GP, not_equal_to )<double> BinaryOpNotEqualTo;

/// less
extern CC_NS( ARM_GP, less )<double> BinaryOpLess;

/// greater
extern CC_NS( ARM_GP, greater )<double> BinaryOpGreater;

/// less_equal
extern CC_NS( ARM_GP, less_equal )<double> BinaryOpLessEqual;

/// greater_equal
extern CC_NS( ARM_GP, greater_equal )<double> BinaryOpGreaterEqual;

/// logical_and
extern CC_NS( ARM_GP, logical_and )<double> BinaryOpLogicalAnd;

/// logical_or
extern CC_NS( ARM_GP, logical_or )<double> BinaryOpLogicalOr;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

