/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file checkinputs.cpp
 *
 *  \brief file to check inputs
 *
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date November 2003
 */

#include "gpbase/checkinputs.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpvector.h"
#include "expt.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Routine: CheckNbSmaller
///	Action : checks that nb1 smaller than nb2 
////////////////////////////////////////////////////

void CheckNbSmaller( double nb1, double nb2, const string& nb1Name, 
        const string& nb2Name, const string& funcName,
        long LINE, char* FILE)
{
	if( nb1 > nb2 )
	{
		CC_Ostringstream os;
		os << "Trying to use " << funcName << " with a shorter " << nb1Name 
			<< " than " << nb2Name << " " << ARM_USERNAME << ": please advise";
		throw Exception(LINE, FILE, ERR_INVALID_ARGUMENT, os.str() );
	}
}


////////////////////////////////////////////////////
///	Routine: CheckPositiveNb
///	Action : Checks that a number is positive
////////////////////////////////////////////////////

void CheckPositiveNb( double nb, const string& nbName, const string& funcName ,
        long LINE, char* FILE)
{
	if( nb < 0 )
	{
		CC_Ostringstream os;
		os << "Trying to use a " << funcName << " with a non strictly positive " << nbName 
			<< " " << ARM_USERNAME  << ": please advise";
		throw Exception(LINE,  FILE, ERR_INVALID_ARGUMENT, os.str() );
	}
}



////////////////////////////////////////////////////
///	Routine: CheckVectorPositiveNb
///	Action : Checks that numbers are positive within the vector!
////////////////////////////////////////////////////

void CheckVectorPositiveNb( const std::vector<double>& vec, 
        const string& vecName, const string& funcName,
        long LINE, char* FILE)
{
	for(int i=0; i<vec.size(); ++i )
	{
		if( vec[i] < 0 )
		{
			CC_Ostringstream os;
			os << "Trying to use " << funcName << " with a non strictly positive " << vecName << ": "
				<< vecName << "[" << i << "] = " << vec[i] << " "
				<< ARM_USERNAME  << ": please advise";
			throw Exception(LINE,  FILE, ERR_INVALID_ARGUMENT, os.str() );
		}
	}

}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

