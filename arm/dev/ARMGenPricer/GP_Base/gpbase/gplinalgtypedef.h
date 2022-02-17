/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gplinalgtypedef.h
 *
 *  \brief typedef for the gp linalg objects
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPINFRA_GPLINALGTYPEDEF_H
#define _INGPINFRA_GPLINALGTYPEDEF_H

#include "port.h"
#include "env.h"
#include <deque>

CC_BEGIN_NAMESPACE( ARM )

template <typename T> class ARM_GP_T_Vector;
template <typename T> class ARM_GP_T_Matrix;
template <typename T> class ARM_GP_T_TriangularMatrix;
template <typename T> class ARM_GP_T_Tensor;

/// double part
typedef ARM_GP_T_Vector<double> ARM_GP_Vector;
typedef ARM_GP_T_Matrix<double> ARM_GP_Matrix;
typedef ARM_GP_T_Tensor<double> ARM_GP_Tensor;
typedef ARM_GP_T_TriangularMatrix<double> ARM_GP_TriangularMatrix;;

/// int part
typedef ARM_GP_T_Vector< int >					        ARM_IntVector;
typedef ARM_GP_T_Vector< string >					    ARM_GP_StrVector;
typedef ARM_GP_T_Matrix< int >					        ARM_IntMatrix;
typedef ARM_GP_T_Tensor< int >					        ARM_IntTensor;
typedef ARM_GP_T_Vector< ARM_GP_StrVector >			    ARM_GP_StrVectorVector;

/// bool part
// FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers
typedef std::deque< bool >							    ARM_BoolVector;
typedef ARM_GP_T_Matrix< bool >                         ARM_BoolMatrix;
typedef ARM_GP_T_Tensor< bool >					        ARM_BoolTensor;


CC_END_NAMESPACE()


#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
