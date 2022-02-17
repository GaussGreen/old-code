/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gplinalgconvert.h
 *
 *  \brief conversion routine from the gp linalg object to 
 *			kernel type linalg object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPINFRA_GPLINALGCONVERT_H
#define _INGPINFRA_GPLINALGCONVERT_H

#include "port.h"
#include "env.h"
#include "gplinalgtypedef.h"
#include "typedef.h"

#include <vector>
CC_USING_NS(std,vector)
#include <string>
CC_USING_NS(std,string)

/// forward declaration
class ARM_Vector;
class ARM_Matrix;

CC_BEGIN_NAMESPACE( ARM )

/// conversion between ARM_Vector and std::vector<double>
extern ARM_Vector		To_ARM_Vector( const std::vector<double>* vec );
extern std::vector<double>	To_ARM_GP_Vector( const ARM_Vector& vec );
extern std::vector<double>	To_ARM_GP_Vector( const ARM_Vector* vec );
extern ARM_Vector*		To_pARM_Vector( std::vector<double>* vec );
extern std::vector<double>*	To_pARM_GP_Vector( ARM_Vector* vec );

/// convertion betwenn ARM_Matrix and ARM_GP_Matrix
extern ARM_Matrix		To_ARM_Matrix( const ARM_GP_Matrix& mat );
extern ARM_GP_Matrix	To_ARM_GP_Matrix( const ARM_Matrix& mat );
extern ARM_Matrix*		To_pARM_Matrix( ARM_GP_Matrix* mat );
extern ARM_GP_Matrix*	To_pARM_GP_Matrix( ARM_Matrix* mat );

/// conversion from std::vector to std::vector<double>
extern std::vector<double>* CreateARMGPVectorFromVECTOR( const vector<double>& vec );
extern std::vector<double>* CreateARMGPMAtuVectorFromStrVECTOR(const vector<string>& vec);
extern std::vector<double>* CreateARMVectorFromXLDATEVECTOR( const vector<double>& vec );

/// conversion between matrices types
extern ARM_TriangularMatrixVector ConvertToTriangularMatrixVector( const ARM_MatrixVector& vec );
extern ARM_MatrixVector ConvertToMatrixVector( const ARM_TriangularMatrixVector& vec );
extern ARM_TriangularMatrixVector CloneAndConvertToTriangularMatrixVector( const ARM_MatrixVector& vec );
extern ARM_MatrixVector CloneAndConvertToMatrixVector( const ARM_TriangularMatrixVector& vec );

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
