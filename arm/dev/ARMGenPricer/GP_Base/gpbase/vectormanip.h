/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: vectormanip.h,v $
 * Revision 1.1  2004/02/06 16:45:06  ebenhamou
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file vectormanip.h
 *
 *  \brief files to handle string argument ... similar
 *			to the interface interglob.cpp file
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_VECTORMANIP_H
#define _INGPBASE_VECTORMANIP_H

#include "port.h"
#include <string>
#include "gplinalgtypedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// merge the two vector assuming they are sorted!
ARM_GP_Vector* MergeSortedVectorNoDuplicates( const ARM_GP_Vector& vec1, const ARM_GP_Vector& vec2 );
ARM_GP_Vector* MergeSortedVectorWithDuplicates( const ARM_GP_Vector& vec1, const ARM_GP_Vector& vec2 );

/// sort using the efficient STL sort algorithm
ARM_GP_Vector* SortSTLBased( const ARM_GP_Vector& vec );

/// remove duplicates
ARM_GP_Vector* VectorUnique(const ARM_GP_Vector& vec );

CC_NS(std,string) VectorToString( const ARM_GP_Vector& vec );

bool ExistsInVector( const ARM_GP_Vector& dateVector, double date );

size_t IdxFromValue( const ARM_GP_Vector& vec, double value);

size_t IdxFromValueWithTol( const ARM_GP_Vector& vec, double value, double tol , bool NoThrow = false);


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
