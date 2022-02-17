/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file checkinputs.h
 *
 *  \brief various function to check inputs
 *
 *	\author  E. Benhamou J.M Prie
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPBASE_CHECKINPUTS_H
#define _INGPBASE_CHECKINPUTS_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "env.h"
#include "port.h"
#include "gplinalgtypedef.h"
#include "ostringstream.h"
#include <string>
CC_USING_NS( std, string )

/// kernel
#include "expt.h"

CC_BEGIN_NAMESPACE( ARM )

/// template version
template< typename T> void CheckVectorIncreasing( const T& vec, const string& vecName, const string& funcName, 
	long LINE = __LINE__, char* FILE = __FILE__  )
{
	for(size_t i=1; i<vec.size(); ++i)
	{
		if( vec[i]< vec[i-1] )
		{
			CC_Ostringstream os;
			os << "Trying to use " << funcName << " with " << vecName << " that are non strictly increasing: " 
				<< vecName << "[" << i-1 << "]: " << vec[i-1] << " > "
				<< vecName << "[" <<  i << "]: " << vec[i]   << " "
				<< ARM_USERNAME  << ": please advise";
			throw Exception(LINE, FILE, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
}


template< typename T> void CheckVectorStrictlyIncreasing( const T& vec, const string& vecName, const string& funcName,
	long LINE = __LINE__, char* FILE = __FILE__  )
{
	for(size_t i=1; i<vec.size(); ++i)
	{
		if( vec[i] <= vec[i-1])
		{
			CC_Ostringstream os;
			os << "Trying to use " << funcName << " with " << vecName << " that are non strictly increasing: " 
				<< vecName << "[" << i-1 << "]: " << vec[i-1] << " > "
				<< vecName << "[" <<  i << "]: " << vec[i]   << " "
				<< ARM_USERNAME  << ": please advise";
			throw Exception(LINE, FILE, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
}

/// validation functions
void CheckNbSmaller( double nb1, double nb2, const string& nb1Name, const string& nb2Name, const string& funcName,
        long LINE = __LINE__, char* FILE = __FILE__ );
void CheckPositiveNb( double nb, const string& nbName, const string& funcName,
        long LINE = __LINE__, char* FILE = __FILE__   );
void CheckVectorPositiveNb( const std::vector<double>* vec, const string& vecName, const string& funcName,
        long LINE = __LINE__, char* FILE = __FILE__  );



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

