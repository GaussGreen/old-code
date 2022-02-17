/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mktdataretriever.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file mktdataretriever.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPCALCULATORS_MKTDATARETRIEVER_H
#define _INGPCALCULATORS_MKTDATARETRIEVER_H

#include "gpbase/removeidentifiedwarning.h"
#include "mktdatacst.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/ostringstream.h"

/// STL
#include <string>
CC_USING_NS(std,string)
#include <vector>
CC_USING_NS(std,vector)

/// ARM Kernel
#include <glob/dates.h>
#include <glob/expt.h>

/// forward declaration in global namespace
class ARM_ZeroCurve;
class ARM_VolLInterpol;
class ARM_VolCube;
class ARM_Model;

CC_BEGIN_NAMESPACE( ARM )




/////////////////////////////////////////////////////////////////
/// \struct ARM_MarketData_Retriever
/// this class is responsible for retrieving data from Summit
/////////////////////////////////////////////////////////////////

/// forward declaration
template <typename T> class ARM_SingletonHolder;

struct ARM_MarketData_RetrieverImp
{
private:
	/// validation of data
	inline void ValidateAsOfDate( const ARM_Date& asOfDate ) const;

	/// conversion function
	inline ARM_Date ConvertDateStringToDate( const string& dateString ) const;
	
	ARM_Date ConvertDateOrMatuString( const string& matuOrDate,
		const string& ccy,
		const ARM_Date& asOfDate,
		int Rule = 0 ) const;
	
	ARM_GP_Vector* ConvertVectorFromMatuOrDateToYearFrac( const vector<string>& vecData, 
		const ARM_Date& asOfDate,
		const string& ccy ) const;

	ARM_GP_Vector* ConvertVectorFromMatuOrDateWRefDateToYearFrac( const vector<string>& vecData, 
		const ARM_Date& asOfDate,
		const ARM_Date& refDate,
		const string& ccy ) const;

	/// various methods to get data
	FILE* GetSummitFile(
		const string& commonFolderName,
		const string& commonFileName,
		const ARM_Date& asOfDate ) const;

	/// function for volatilities Warning the ccy and the asof are required to create the object!
	ARM_VolLInterpol* GetVolFromSummit_Common(
		const string& commonFolderName,
		const string& commonFileName,
		const string& ccy,
		const ARM_Date& asOfDate,
		int volType ) const;

	ARM_VolLInterpol* GetVolSmileFromSummit(
		const string& index,
		const string& ccy,
		const string& cvName,
		const ARM_Date& asOfDate,
		const string& volMktType,
		const string& matuIndex ) const;

	ARM_VolLInterpol* GetVolATMFromSummit(
		const string& index,
		const string& ccy,
		const string& cvName,
		const ARM_Date& asOfDate,
		const string& volMktType ) const;


	ARM_VolCube* GetVolCubeFromSummit(
		const string& index,
		const string& ccy,
		const string& cvName,
		const ARM_Date& asOfDate,
		const string& volMktType ) const;

	/// to forbid client from using it except for the singleton holder
	ARM_MarketData_RetrieverImp() {};
	~ARM_MarketData_RetrieverImp() {};
	/// and non implemented copy and assignment operator to avoid duplication
	ARM_MarketData_RetrieverImp( const ARM_MarketData_RetrieverImp& rhs );
	ARM_MarketData_RetrieverImp& operator=( const ARM_MarketData_RetrieverImp& rhs );

	friend class ARM_SingletonHolder<ARM_MarketData_RetrieverImp>;

public:
	/// functions to get data
	ARM_ZeroCurve* GetZCCurveFromSummit( 
		const string& index,
		const string& ccy,
		const string cvName,
		const ARM_Date& asOfDate ) const;

	ARM_VolLInterpol* GetVolFromSummit(
		const string& index,
		const string& ccy,
		const string& cvName,
		const ARM_Date& asOfDate,
		ARM_MarketData::VolMktType volMktType,
		ARM_MarketData::VolType volType ) const;

	ARM_Model* GetVolMktModel( 
		const string& indexName, 
		const string& ccy, 
		const string& cvName, 
		const ARM_Date& asOf, 
		const string& volMktType,
		const string& source ) const;
};

extern ARM_SingletonHolder<ARM_MarketData_RetrieverImp> ARM_TheMarketDataRetriever;

//////////////////////////////////////////////////////////
/// inline code part for fast access
//////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
///	Action : checks that the asOfDate is correctly set up
/////////////////////////////////////////////////////////////////
inline void ARM_MarketData_RetrieverImp::ValidateAsOfDate(	const ARM_Date& asOfDate ) const
{
	/// some validation:
	/// 1) are we in the future?
	static ARM_Date Today;
	
	/// test that we are not after today!
	if( asOfDate > Today )
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << ": trying to access a curve with asOfDate " << asOfDate.toString()
			<< " while today= " << Today.toString();
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
}


/////////////////////////////////////////////////////////////////
///	Action : convert a dateString according to the default setting
/////////////////////////////////////////////////////////////////
inline ARM_Date ARM_MarketData_RetrieverImp::ConvertDateStringToDate( const string& dateString ) const
{
	if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
		return ARM_Date( (char*) dateString.c_str(),"MM/DD/YYYY");
	else
		return ARM_Date( (char*) dateString.c_str(),"DD/MM/YYYY");
}



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
