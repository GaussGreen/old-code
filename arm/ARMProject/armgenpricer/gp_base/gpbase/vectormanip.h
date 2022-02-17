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
std::vector<double>* MergeSortedVectorNoDuplicates( const std::vector<double>* vec1, const std::vector<double>* vec2 );
std::vector<double>* MergeSortedVectorWithDuplicates( const std::vector<double>* vec1, const std::vector<double>* vec2 );

/// sort using the efficient STL sort algorithm
std::vector<double>* SortSTLBased( const std::vector<double>* vec );

/// remove duplicates
std::vector<double>* VectorUnique(const std::vector<double>* vec );

CC_NS(std,string) VectorToString( const std::vector<double>* vec );

bool ExistsInVector( const std::vector<double>* dateVector, double date );

size_t IdxFromValue( const std::vector<double>* vec, double value);

size_t IdxFromValueWithTol( const std::vector<double>* vec, double value, double tol , bool NoThrow = false);


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
