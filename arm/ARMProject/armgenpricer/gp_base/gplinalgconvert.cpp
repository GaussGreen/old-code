/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gplinalgconvert.cpp
 *
 *  \brief file to convert linalg object to other linalg object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpbase/gplinalgconvert.h"

#include "gpbase/gpmatrix.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/gpvector.h"
#include "gpbase/stringconvert.h"
#include "gpbase/datemanip.h"

#include "linalg.h"


CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////
/// conversion routines from ARM_Vector to std::vector<double> and vice versa
/////////////////////////////////////////////////////////

ARM_Vector  To_ARM_Vector( const std::vector<double>& vec )
{
	return ARM_Vector( vec ); 
}

std::vector<double> To_ARM_GP_Vector( const ARM_Vector& vec )
{ 
	return std::vector<double>( vec.begin(), vec.end() );
}

std::vector<double> To_ARM_GP_Vector( const ARM_Vector* vec )
{
	if (!vec)
		return std::vector<double>();
	else
		return std::vector<double>( vec->begin(), vec->end());
}

/// pointor version (beware of memory leak as this creates a new copy!)
ARM_Vector*  To_pARM_Vector( std::vector<double>* vec )
{ 
	return vec? new ARM_Vector(vec->size(),(double*)(&(*vec->begin()))): NULL;
}

std::vector<double>* To_pARM_GP_Vector( ARM_Vector* vec )
{ 
	return vec? new std::vector<double>(vec->begin(),vec->end()): NULL; 
}



/////////////////////////////////////////////////////////
/// conversion routines from ARM_Matrix to ARM_GP_Matrix and vice versa
/////////////////////////////////////////////////////////

ARM_Matrix  To_ARM_Matrix( const ARM_GP_Matrix& mat )
{
// FIXMEFRED: mig.vc8 (22/05/2007 15:57:58):explicit cast
	return ARM_Matrix( mat.rows(), mat.cols(), (double*) &(*mat.begin()) ); 
}

ARM_GP_Matrix To_ARM_GP_Matrix( const ARM_Matrix& mat )
{ 
	return ARM_GP_Matrix( (size_t) mat.GetNumLines(), (size_t) mat.GetNumCols(), (double*) &mat.Elt(0,0) );
}

/// pointor version (beware of memory leak as this creates a new copy!)
ARM_Matrix*  To_pARM_Matrix( ARM_GP_Matrix* mat )
{ 
// FIXMEFRED: mig.vc8 (22/05/2007 15:57:58):explicit cast
	return mat? new ARM_Matrix(mat->rows(), mat->cols(), (double*)&(*mat->begin())): NULL;
}

ARM_GP_Matrix* To_pARM_GP_Matrix( ARM_Matrix* mat )
{ 
	return mat? new ARM_GP_Matrix( (size_t) mat->GetNumLines(), (size_t) mat->GetNumCols(), (double*) &mat->Elt(0,0) ): NULL; 
}


/////////////////////////////////////////////////////////
/// conversion from std::vector to std::vector<double>
/////////////////////////////////////////////////////////
std::vector<double>* CreateARMGPVectorFromVECTOR( const vector<double>& vec )
{
	return new std::vector<double>( vec );
}

std::vector<double>* CreateARMGPMAtuVectorFromStrVECTOR(const vector<string>& vec)
{
	if(vec.empty()) 
		return NULL;
	else
	{
		size_t vecSize = vec.size();
		std::vector<double>* res = new std::vector<double>(vecSize);
		for( size_t i=0; i<vecSize; ++i )
			(*res)[i] = StringMaturityToYearTerm( vec[i] );
		return res;
	}
}

std::vector<double>* CreateARMVectorFromXLDATEVECTOR( const vector<double>& vec )
{
	if(vec.empty()) 
		return NULL;
	else
	{
		size_t vecSize = vec.size();
		std::vector<double>* res = new std::vector<double>(vecSize);
		for( size_t i=0; i<vecSize; ++i )
			(*res)[i] = ConvertXLDateToJulian( vec[i] );
		return res;
	}	
}



/////////////////////////////////////////////////////////
/// matrix conversion routines
/////////////////////////////////////////////////////////
ARM_TriangularMatrixVector ConvertToTriangularMatrixVector( const ARM_MatrixVector& vec )
{
	ARM_TriangularMatrixVector result( vec.size() );
	for( size_t i=0; i<vec.size(); ++i )
		result[i] = (ARM_GP_TriangularMatrix*) vec[i];
	return result;
}


ARM_TriangularMatrixVector CloneAndConvertToTriangularMatrixVector( const ARM_MatrixVector& vec )
{
	ARM_TriangularMatrixVector result( vec.size() );
	for( size_t i=0; i<vec.size(); ++i )
		result[i] = new ARM_GP_TriangularMatrix( *vec[i] );
	return result;
}


ARM_MatrixVector ConvertToMatrixVector( const ARM_TriangularMatrixVector& vec )
{
	ARM_MatrixVector result( vec.size() );
	for( size_t i=0; i<vec.size(); ++i )
		result[i] = (ARM_GP_Matrix*) vec[i];
	return result;
}

ARM_MatrixVector CloneAndConvertToMatrixVector( const ARM_TriangularMatrixVector& vec )
{
	ARM_MatrixVector result( vec.size() );
	for( size_t i=0; i<vec.size(); ++i )
	{
		result[i] = new ARM_GP_Matrix( vec[i]->rows(), vec[i]->cols() );

		for( size_t j=0; j<vec[i]->rows(); ++j )
			for( size_t k=0; k<vec[i]->cols(); ++k )
				result[i]->Elt(i,k)= vec[i]->Elt(j,k);
	}
	return result;
}




CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
