/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gptrigomatrix.h
 *
 *  \brief a matrix with trigonometric coefficient used
 *  for correlation
 *	\author  Richard GUILLEMOT
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPINFRA_GPTRIGOMATRIX_H
#define _INGPINFRA_GPTRIGOMATRIX_H

#include "port.h"
#include "env.h"
#include "gpmatrix.h"
#include "gplinalgtypedef.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_GP_Matrix* TrigoMatrix(size_t n, double alpha);

CC_END_NAMESPACE()

#endif