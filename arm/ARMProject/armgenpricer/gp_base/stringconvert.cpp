/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: stringconvert.cpp,v $
 * Revision 1.1  2003/10/08 16:39:43  ebenhamou
 * Initial revision
 *
 */

/*! \file stringconvert.cpp
 *
 *  \brief file to convert string to double and vice versa!
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/stringconvert.h"
#include "gpbase/ostringstream.h"
#include <cstdlib>	/// for toupper
#include <cmath>	/// for fabs
#include "expt.h"
#include "armglob.h"
//#include <util/fromto.h>
//#include <inst/irindex.h>

CC_BEGIN_NAMESPACE( ARM )

/// function to convert a string representing a maturity to a year fraction
double StringMaturityToYearTerm( const string& maturity )
{
	
    int nb;
    char unit;
	sscanf(maturity.c_str(), "%d%c", &nb, &unit);
	unit = toupper(unit);

    switch (unit)
    {
        case 'M' : return(double(nb)/12.0); //*0.0833333333333333;			/// = 1.0/12.0;
        case 'Y' : return double(nb);								/// = 1.0;
        case 'D' : case 'B' : return (double(nb)*0.00273972602739726);	/// = 1.0/365.0;
        case 'W' : return (double(nb)*0.0191780821917808);			/// = 7.0/365.0
        default  : 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"invalid maturity string should be like nb[d/w/m/y]");
    }
}

//////////////////////////////////////////
/// returns true if manage to convert
/// else return false
//////////////////////////////////////////
bool YearTermToStringMaturityHelper( double d, string& res, double period, const string& periodString )
{
	const double precision = 1e-1;

	CC_Ostringstream os;

	/// priority to year then to month then week then to day to day
	if( fabs( period * d - int( period * d ) ) < precision )
	{
		os << int( period  *d ) << periodString;
		res = os.str();
		return true;
	}
	else
		return false;
}


//////////////////////////////////////////
/// because of potential ambiguity 
/// priority to year then to month then week then to day to day
/// choice is that for more than 2 year
/// we round to the corresponding year
/////////////////////////////////////////
string YearTermToStringMaturity( double d )
{
	string result;
	CC_Ostringstream os;

	if( d > 2.0 )
	{
		os << d << "Y";
		return os.str();
	}

	if( YearTermToStringMaturityHelper( d, result, 1.0, "Y" ) )
		return result;
	if( YearTermToStringMaturityHelper( d, result, 12.0, "M" ) )
		return result;
	if( YearTermToStringMaturityHelper( d, result, 52.1428571428571, "W" ) ) /// = 365.0/52.0
		return result;
	if( YearTermToStringMaturityHelper( d, result, 365.0, "D" ) )
		return result;
	/// if did not manage to convert throw an exception
	throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "Invalid double period");
}



/// from the term and the index as string, 
/// return the maturity type as an int define (of type ARM_INDEX_TYPE)
int FromIndexAndTermToIndexType( const string& term, const string& index )
{
    if (index == "FIXED")
        return K_FIXED;

    if (index == "EUR")
    {
        if (term == "1M")
            return K_EURIBOR1M;
        if (term == "3M")
            return K_EURIBOR3M;
        if (term == "6M")
            return K_EURIBOR6M;
        if (term == "1Y")
            return K_EURIBOR1Y;
		if (term == "12M")
            return K_EURIBOR1Y;

        throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
            "Not a good term for Euribor");
    }
    else 
		if (index == "PIBOR")
		{
			if (term =="1M")
				return K_PIBOR1M;
			if (term == "3M")
				return K_PIBOR3M;
			if (term =="6M")
				return K_PIBOR6M;
			if (term == "1Y")
				return K_PIBOR1Y;
			throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
				"Not a good term for Pibor");
		}
		else
		{
			/// For the other index, use libor constant
			if (term =="1M")
				return K_LIBOR1M;
			if (term == "3M")
				return K_LIBOR3M;
			if (term == "6M")
				return K_LIBOR6M;
			if (term == "1Y")
				return K_LIBOR1Y;
			/// other cases are errors!
			throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
				"Not a good term for Libor");
		}
}


/// from the term and the index as string, 
/// return the maturity type as an int define (of type ARM_INDEX_TYPE)
string FromIndexTypeToTermAndType( int indexType, int& type )
{
	CC_Ostringstream tenorDesc;

	//// Fixed
	//if (IsFixedIndex((ARM_INDEX_TYPE)indexType))
	//{
	//	type = K_FIXED;
	//}
	//// Libor
	//else if (IsLiborIndex((ARM_INDEX_TYPE)indexType))
	//{
	//	type = K_LIBOR;
	//	ARM_IRIndex index((ARM_INDEX_TYPE)indexType);
	//	double tenor = index.GetYearTerm();
	//	int nmonth = int(tenor * 12);
	//	tenorDesc  << nmonth  << "m" ;
	//}
	///// CMS case
	//else if (IsCMSIndex((ARM_INDEX_TYPE)indexType))
	//{
	//	type = K_CMS;
	//	tenorDesc  << indexType  - K_CMS1 + 1  << "y" ;
	//}
	//else
	//{
	//	throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
	//			"The index type is not valid.");
	//}

	return tenorDesc.str();

}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

