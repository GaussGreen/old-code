/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curve.cpp
 *  \brief file for the definition of templated curves
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date May 2004
 */

#include "gpbase/curve.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
/// code part
////////////////////////////////////////////////////
/// function to paliate the pb of ostringstream not working in debug mode
int Sprintf( char* msg, double value )
{	return sprintf( msg, "%f", value ); }


/// function to paliate the pb of ostringstream not working in debug mode
int Sprintf( char* msg, const vector<double>& vec )
{
	int pos = 0;
	for( size_t i=0; i<vec.size(); ++i )
	{	pos += sprintf( msg+pos, "%f\t", vec[i] ); }
	return pos;
}


CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
